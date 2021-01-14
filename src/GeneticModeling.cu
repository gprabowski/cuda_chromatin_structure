#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_fp16.h>

#include <thrust/fill.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>

#include <time.h>
#include <iostream>

#include <stdio.h>
#include "../include/LooperSolver.h"

#define ITERATIONS 1000
#define ISLAND_SIZE 32
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
#define cucheck_dev(call)                                   \
{                                                           \
  cudaError_t cucheck_err = (call);                         \
  if(cucheck_err != cudaSuccess) {                          \
    const char *err_str = cudaGetErrorString(cucheck_err);  \
    printf("%s (%d): %s\n", __FILE__, __LINE__, err_str);   \
    assert(0);                                              \
  }                                                         \
}

typedef struct half3 {
	__half x;
	__half y;
	__half z;
} half3;

typedef struct gpu_settings {
	float densityInfluence;
	float densityScale;
	float densityWeight;
	double maxTempHeatmapDensity;
	float dtTempHeatmapDensity;
	float dtTempHeatmap;
	double eps;
	float MCstopConditionImprovementHeatmapDensity;
	float MCstopConditionImprovementHeatmap;
	bool use2D;
	float tempJumpCoefHeatmapDensity;
	float tempJumpCoefHeatmap;
	float tempJumpScaleHeatmapDensity;
	float tempJumpScaleHeatmap;
	float maxTempHeatmap;
} gpu_settings;

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
	if (code != cudaSuccess) {
    	fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      	if (abort) exit(code);
   	}
}

__global__ void setupKernel(curandState * state, time_t seed) {
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    curand_init(seed, tid, 0, &state[tid]);
}

__device__ void bitonicSort(__half *dev_values, short *indices, int islandSize) {
	short i, ixj; /* Sorting partners: i and ixj */

	i = threadIdx.x;
	indices[i] = i;
	
	for (short k = 2; k <= islandSize; k <<= 1) {
		for (short j = k >> 1; j > 0; j >>= 1) {
			ixj = i ^ j;

			/* The threads with the lowest ids sort the array. */
			if (ixj > i) {
				if ((i&k) == 0) {
					/* Sort ascending */
					if (__hgt(dev_values[i], dev_values[ixj])) {
						/* exchange(i,ixj); */
						__half temp = dev_values[i];
						dev_values[i] = dev_values[ixj];
						dev_values[ixj] = temp;
						indices[i] = ixj;
						indices[ixj] = i;
					}
				}
				if ((i&k) != 0) {
					/* Sort descending */
					if (__hgt(dev_values[ixj], dev_values[i])) {
						/* exchange(i,ixj); */
						__half temp = dev_values[i];
						dev_values[i] = dev_values[ixj];
						dev_values[ixj] = temp;
						indices[i] = ixj;
						indices[ixj] = i;
					}
				}
			}
			__syncwarp();
		}
	}
}

__device__ void subtractFromVector3(half3 *destination, half3 *value1, half3 *value2) {
	destination->x = __hsub(value1->x, value2->x);
	destination->y = __hsub(value1->y, value2->y);
	destination->z = __hsub(value1->z, value2->z);
}

__device__ __half magnitude(half3 *vector) {
	__half sum = 0;
	
	sum = __hadd(sum, __hmul(vector->x, vector->x));
	sum = __hadd(sum, __hmul(vector->y, vector->y));
	sum = __hadd(sum, __hmul(vector->z, vector->z));

	// return sqrtf(__half2float(sum));
	return hsqrt(sum);
}

__device__ void getChromosomeHeatmapBoundary(int p, int &start, int &end, int* gpu_heatmap_chromosome_boundaries, int chromosomeBoundariesSize) {
	if (chromosomeBoundariesSize == 0) return;
	if (p < 0 || p > gpu_heatmap_chromosome_boundaries[chromosomeBoundariesSize - 1]) return;

	for (size_t i = 0; i+1 < chromosomeBoundariesSize; ++i) {
		if (gpu_heatmap_chromosome_boundaries[i] <= p && p < gpu_heatmap_chromosome_boundaries[i+1]) {
			start = gpu_heatmap_chromosome_boundaries[i];
			end = gpu_heatmap_chromosome_boundaries[i+1]-1;
			return;
		}
	}
}

__device__ float calcScoreHeatmapSingleActiveRegion(
    int moved,
    half3 * clusters_positions,
    int* gpu_heatmap_chromosome_boundaries,
    float * heatmap_dist,
    int heatmapSize,
    int heatmapDiagonalSize,
    int activeRegionSize,
    int chromosomeBoundariesSize
) {
    float err = 0.0, cerr;
	half3 temp;
	float d;
    size_t n = activeRegionSize;
    if (heatmapSize != n) { printf("heatmap sizes mismatch, dist size=%d, active region=%zu", heatmapSize, n); return 0.0; }
    int st = 0;
    int end = n - 1;
    
    if(chromosomeBoundariesSize > 0) {
        getChromosomeHeatmapBoundary(moved, st, end, gpu_heatmap_chromosome_boundaries, chromosomeBoundariesSize);
    }

    for (int i = st; i <= end; ++i) {
        if (abs(i-moved) >= heatmapDiagonalSize) {
            if (heatmap_dist[i * heatmapSize + moved] < 1e-6) continue;	// ignore 0 values
            half3* temp_one = clusters_positions + i;
            half3* temp_two = clusters_positions + moved;
            subtractFromVector3(&temp, temp_one, temp_two);
            d = magnitude(&temp);
            cerr = (d - heatmap_dist[i * heatmapSize + moved]) / heatmap_dist[i * heatmapSize + moved];
            err += cerr * cerr;
        }
    }
    return err;
}

__global__ void calculateScoreHeatmapActiveRegion(
	half3 * clusters_positions,
	int* gpu_heatmap_chromosome_boundaries, 
	float * heatmap_dist,
	__half * scores,
    int heatmapSize,
    int heatmapDiagonalSize,
    int activeRegionSize,
    int chromosomeBoundariesSize
) {
	short x = threadIdx.x;
	short y = threadIdx.y;
	short z = threadIdx.z;

	// TODO: chromosomeBoundaries (is it even needed?)

	if (x == y || abs(x-y) < heatmapDiagonalSize) return;
	
	__half heatmapValue = __float2half(heatmap_dist[x * heatmapSize + y]);

	if (__hgt(1e-6, heatmapValue)) return; // ignore 0 values

	short xId = blockDim.x * blockIdx.x + threadIdx.x + (z * activeRegionSize);
	short yId = blockDim.y * blockIdx.x + threadIdx.y + (z * activeRegionSize);
	
	half3 temp;
	half3* temp_one = clusters_positions + xId;
    half3* temp_two = clusters_positions + yId;
	subtractFromVector3(&temp, temp_one, temp_two);
	
    __half d = magnitude(&temp);
	__half cerr = __hdiv(__hsub(d, heatmapValue), heatmapValue);
	
	// could be replaced with reduction ( O(n) -> O(log(n)) )
	atomicAdd(&scores[blockIdx.x * ISLAND_SIZE + z], __hmul(cerr, cerr));
}

__device__ void calculateFitness(
	half3 * populationGlobal,
	half3 * populationShared,
	int * heatmap_chromosome_boundaries, 
	float * heatmap_dist,
	__half * scores,
	__half * fitness,
	const int heatmapSize,
	const int heatmapDiagonalSize,
	const int clusterSize,
	const int chromosomeBoundariesSize
) {
	int threadIndex = blockDim.x * blockIdx.x + threadIdx.x;

	dim3 block(clusterSize, clusterSize, ISLAND_SIZE);

	for(int i = 0; i < clusterSize; ++i) {
		populationGlobal[threadIndex * clusterSize + i] = populationShared[threadIdx.x * clusterSize + i];
	}
	
	if(threadIndex == 0) {
		calculateScoreHeatmapActiveRegion<<<gridDim.x, block>>>(
			populationGlobal,
			heatmap_chromosome_boundaries, 
			heatmap_dist,
			scores,
			heatmapSize,
			heatmapDiagonalSize,
			clusterSize,
			chromosomeBoundariesSize
		);
		cucheck_dev( cudaPeekAtLastError() );
		cucheck_dev( cudaDeviceSynchronize() );
	}
	__syncthreads();
	
	fitness[threadIdx.x] = scores[blockIdx.x * ISLAND_SIZE + threadIdx.x];
	__syncthreads();
}

// ACTION ITEMS:
// 1. DONE - use __half to allow for 2x larger population in a block
// 2. try doing crossover on one thread from each pair only (less barriers, but one thread idle)
// 		is it even possible?
// 3. optimize migration by utilizing more than one thread for data transfer
// 4. Since we probably won't be able to fit more than 32 elements per block, we can use warp-level operations
// 5. DONE - OF COURSE WE CAN - Can we do better with the calcScoreHeatmapActiveRegion function?
// 6. We could utilize Tensor Cores for offspring calculations
// 7. Test thrust reduction on the whole population vs local warp reductions => thrust on smaller load
// 8. Test half3 struct vs flat 1D array of __half
// 9. Vectorized operations on half3? __align__(8) - do we have enough SM?
// 10. Use constant or texture memory for read-only arrays
// 11. DONE - merge positions and active region to use less shared memory
// 12. Randomize input data (atm all the solutions are the same)
// 13. Verify if we allocate correct number of elements (population for chr14 has size 1878, but 58 * 32 = 1856)
__global__ void geneticHeatmap(
	curandState * state,
	half3 * populationGlobal,
	int * heatmap_chromosome_boundaries,
	float * heatmap_dist,
	__half * scores,
	gpu_settings settings,
    const int numberOfIslands,
	const int clusterSize,
	const int migrationConstant,
    const int heatmapSize,
    const int heatmapDiagonalSize,
	const int chromosomeBoundariesSize
) {
	int threadIndex = blockDim.x * blockIdx.x + threadIdx.x;
	if(threadIdx.x >= ISLAND_SIZE) return;

	extern __shared__ half3 population[];
	__shared__ __half fitness[ISLAND_SIZE];
	__shared__ short selectedIndices[ISLAND_SIZE];
	__shared__ __half weightsAndCrossover[ISLAND_SIZE];
	__shared__ short sortedIndices[ISLAND_SIZE];

	half3 tempChild[200];
	int crossoverIdx;
	half3 *parent1, *parent2;
	half3 *localChromosome;
	__half a, gauss;
	float mutationProbability;

	curandState localState = state[threadIndex];

	for(int i = 0; i < clusterSize; ++i) {
		population[threadIdx.x * clusterSize + i] = populationGlobal[threadIndex * clusterSize + i];
	}

	int i = 1;

	while(true) {
		localChromosome = &population[threadIdx.x * clusterSize];

		// fitness evaluation
		calculateFitness(populationGlobal, population, heatmap_chromosome_boundaries, heatmap_dist,
			scores, fitness, heatmapSize, heatmapDiagonalSize, clusterSize, chromosomeBoundariesSize
		);

		if(i == ITERATIONS) break;

		// MIGRATION - START
		if(i % 100 == 0) {
			// fitness[threadIdx.x] = __float2half(curand_uniform(&localState));
			bitonicSort(fitness, sortedIndices, ISLAND_SIZE);
			
			if(threadIdx.x == 0) {
				// write M best chromosomes to neighbouring island on the left
				int islandId = blockIdx.x == 0 ? numberOfIslands - 1 : blockIdx.x - 1;

				for(int j = 0; j < migrationConstant; ++j) {
					for(int k = 0; k < clusterSize; ++k) {

						// sortedIndices[j] in global?
						populationGlobal[islandId * ISLAND_SIZE * clusterSize + sortedIndices[j] * clusterSize + k] 
							= population[sortedIndices[j] * clusterSize + k];
					}
				}
				
				// replace M worst island chromosomes with M best chromosomes from neigbouring island to the right
				islandId = blockIdx.x == numberOfIslands - 1 ? 0 : blockIdx.x + 1;

				for(int j = ISLAND_SIZE - migrationConstant; j < ISLAND_SIZE; ++j) {
					for(int k = 0; k < clusterSize; ++k) {
						population[sortedIndices[j] * clusterSize + k] 
							= populationGlobal[islandId * ISLAND_SIZE * clusterSize + sortedIndices[j] * clusterSize + k];
					}
				}
			}
			__syncwarp();
 
			localChromosome = &population[threadIdx.x * clusterSize];

			// fitness evaluation
			calculateFitness(populationGlobal, population, heatmap_chromosome_boundaries, heatmap_dist,
				scores, fitness, heatmapSize, heatmapDiagonalSize, clusterSize, chromosomeBoundariesSize
			);
		}
		// MIGRATION - END

		// selection and crossover
		crossoverIdx = curand(&localState) % ISLAND_SIZE;

		selectedIndices[threadIdx.x] = __hgt(fitness[threadIdx.x], fitness[crossoverIdx]) ? threadIdx.x : crossoverIdx;
		weightsAndCrossover[threadIdx.x] = __float2half(curand_uniform(&localState));
		__syncwarp();

		crossoverIdx = threadIdx.x % 2 == 0 ? threadIdx.x : threadIdx.x - 1;

		if(__hgt(__float2half(0.7), weightsAndCrossover[crossoverIdx])) {
			a = weightsAndCrossover[crossoverIdx + 1];
			crossoverIdx = threadIdx.x % 2 == 0 ? 1 : -1;  // reuse the variable

			parent1 = &population[selectedIndices[threadIdx.x] * clusterSize];
			parent2 = &population[selectedIndices[(int)threadIdx.x + crossoverIdx] * clusterSize];
			__half oneMinusA_fp16 = __hsub(__float2half(1.0), a);

			for(int j = 0; j < clusterSize; ++j) {
				// Offspring1 = a * Parent1 + (1- a) * Parent2
				// Offspring2 = (1 - a) * Parent1 + a * Parent2
				tempChild[j].x = __hadd(__hmul(a, parent1[j].x), __hmul(oneMinusA_fp16, parent2[j].x));
				tempChild[j].y = __hadd(__hmul(a, parent1[j].y), __hmul(oneMinusA_fp16, parent2[j].y));
				tempChild[j].z = __hadd(__hmul(a, parent1[j].z), __hmul(oneMinusA_fp16, parent2[j].z));
			}
		}
		__syncwarp();
		
		// mutation and store back to shared memory
		gauss = __float2half(curand_normal(&localState));
		mutationProbability = curand_uniform(&localState);

		for(int j = 0; j < clusterSize; ++j) {
			localChromosome[j] = tempChild[j];

			if(mutationProbability < 0.05) {
				localChromosome[j].x = __hadd(localChromosome[j].x, gauss);
				localChromosome[j].y = __hadd(localChromosome[j].y, gauss);
				localChromosome[j].z = __hadd(localChromosome[j].z, gauss);
			}
		}
		__syncwarp();
		++i;
	}
}

float LooperSolver::ParallelGeneticHeatmap() {
	const int M = 0.2 * active_region.size();
    const int blocks = 256;  // const int blocks = Settings::numberOfBlocks;
	const int threads = 32;  // const int threads = Settings::numberOfThreads;
    const int numberOfIndividuals = blocks * threads;

	float score, score_density;
    int size = active_region.size();
    int heatmapSize = heatmap_dist.size;

	if (size <= 1) return 0.0;	// there is nothing to do

	gpuErrchk( cudaDeviceSetCacheConfig(cudaFuncCachePreferShared) );

	// calc initial score
	score = calcScoreHeatmapActiveRegion();	// score from heatmap

	// TODO: Do we need this? It's not used anywhere in the algorithm
	score_density = calcScoreDensity();
    
    struct gpu_settings settings;
    settings.tempJumpScaleHeatmap = Settings::tempJumpScaleHeatmap;
    settings.tempJumpCoefHeatmap = Settings::tempJumpCoefHeatmap;
    settings.use2D = Settings::use2D;
	
    thrust::host_vector<half3> h_clusters_positions(size); // just one initial state to copy across
    thrust::host_vector<float> h_scores(numberOfIndividuals);

    for(int i = 0; i < size; ++i) {
        h_clusters_positions[i].x = __float2half(clusters[active_region[i]].pos.x);
        h_clusters_positions[i].y = __float2half(clusters[active_region[i]].pos.y);
        h_clusters_positions[i].z = __float2half(clusters[active_region[i]].pos.z);
    }

    thrust::device_vector<int> d_heatmap_chromosome_boundaries(heatmap_chromosome_boundaries);
    thrust::device_vector<half3> d_clusters_positions(numberOfIndividuals * size);
    thrust::device_vector<__half> d_scores(numberOfIndividuals, 0.0);
    thrust::device_vector<float> d_heatmap_dist(heatmapSize * heatmapSize);

    // TODO: make sure it's copying in the correct order
    for(int i = 0; i < heatmapSize; ++i) {
        thrust::copy(heatmap_dist.v[i], heatmap_dist.v[i] + heatmapSize, d_heatmap_dist.begin() + heatmapSize * i);
    }

	for(int i = 0; i < numberOfIndividuals; ++i) {
        thrust::copy(h_clusters_positions.begin(), h_clusters_positions.end(), d_clusters_positions.begin() + i * size);
    }

	curandState * devStates;
	gpuErrchk(cudaMalloc((void**)&devStates, numberOfIndividuals * sizeof(curandState)));

    // cuRand initialization
	setupKernel<<<blocks, threads>>>(devStates, time(NULL));
	gpuErrchk( cudaDeviceSynchronize() );

	output(2, "initial score: %lf (density=%lf)\n", score, score_density);

	geneticHeatmap<<<blocks, threads, threads * size * sizeof(half3)>>>(
		devStates, //curandState * state,
		thrust::raw_pointer_cast(d_clusters_positions.data()), //__half * populationGlobal,
		thrust::raw_pointer_cast(d_heatmap_chromosome_boundaries.data()), //int * heatmap_chromosome_boundaries,
		thrust::raw_pointer_cast(d_heatmap_dist.data()), //__half * heatmap_dist,
		thrust::raw_pointer_cast(d_scores.data()), //__half * scores,
		settings, //gpu_settings settings,
		blocks, //const int numberOfIslands,
		size, // 3 * active region
		M,
		heatmapSize,
		heatmap_dist.diagonal_size,
		heatmap_chromosome_boundaries.size()
	);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	int resultIndex = thrust::min_element(thrust::device, d_scores.begin(), d_scores.end()) - d_scores.begin();

	auto iterStart = d_clusters_positions.begin() + resultIndex * size;
    auto iterEnd = d_clusters_positions.begin() + resultIndex * size + size;
	
	thrust::copy(iterStart, iterEnd, h_clusters_positions.begin());

	h_scores = d_scores;
	score = h_scores[resultIndex];

    for(int i = 0; i < size; ++i) {
        clusters[active_region[i]].pos.x = h_clusters_positions[i].x;
        clusters[active_region[i]].pos.y = h_clusters_positions[i].y;
        clusters[active_region[i]].pos.z = h_clusters_positions[i].z;
	}
	
	std::cout << "Gentic Algorithm score: " << score << std::endl;

	cudaFree(devStates);
	return score;
}