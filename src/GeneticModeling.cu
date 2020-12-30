#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_fp16.h>

#include <time.h>

#include "../include/LooperSolver.h"

#define ISLAND_SIZE 32
#define ITERATIONS 1000
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

typedef __half FP16;

typedef struct __align__(6) {
	FP16 x;
	FP16 y;
	FP16 z;
 } half3;

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ void bitonicSort(FP16 *dev_values, ushort *indices, int islandSize) {
	unsigned int i, ixj; /* Sorting partners: i and ixj */
	int j, k; 

  	i = threadIdx.x;
	ixj = i ^ j;
	
	for (k = 2; k <= islandSize; k <<= 1) {
		for (j = k >> 1; j > 0; j >>= 1) {
			ixj = i ^ j;

			/* The threads with the lowest ids sort the array. */
			if ((ixj) > i) {
				if ((i&k) == 0) {
					/* Sort ascending */
					if (__hgt(dev_values[i], dev_values[ixj])) {
						/* exchange(i,ixj); */
						FP16 temp = dev_values[i];
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
						FP16 temp = dev_values[i];
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

__device__ FP16 calcScoreHeatmapActiveRegion(
    int moved, 
    FP16* clusters_positions, 
    int* active_region, 
    int* gpu_heatmap_chromosome_boundaries, 
    FP16 * heatmap_dist,
    int heatmapSize,
    int heatmapDiagonalSize,
    int activeRegionSize,
    int chromosomeBoundariesSize
) {
    float err = 0.0;
    size_t n = activeRegionSize;
	if (moved == -1) {
        for (size_t i = 0; i < n; ++i) {
            err += calcScoreHeatmapSingleActiveRegion(i, clusters_positions, active_region, gpu_heatmap_chromosome_boundaries, 
                heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize);
        }
	}
	else {
        err = calcScoreHeatmapSingleActiveRegion(moved, clusters_positions, active_region, gpu_heatmap_chromosome_boundaries, 
                heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize);
	}
	return __float2half(err);
}

// ACTION ITEMS:
// 1. DONE - use FP16 to allow for 2x larger population in a block
// 2. try doing crossover on one thread from each pair only (less barriers, but one thread idle)
// 		is it even possible?
// 3. optimize migration by utilizing more than one thread for data transfer
// 4. Since we probably won't be able to fit more than 32 elements per block, we can use warp-level operations
// 5. Can we do better with the calcScoreHeatmapActiveRegion function?
// 6. We could utilize Tensor Cores for offspring calculations
// 7. Test thrust reduction on the whole population vs local warp reductions => thrust on smaller load
// 8. Test half3 struct vs flat 1D array of __half
__global__ void geneticHeatmap(
	curandState * state,
	FP16 * populationGlobal,
	int * heatmap_chromosome_boundaries,
	FP16 * heatmap_dist,
	FP16 * scores,
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

	extern __shared__ FP16 population[];
	__shared__ FP16 fitness[ISLAND_SIZE];
	__shared__ ushort selectedIndices[ISLAND_SIZE];
	__shared__ FP16 weightsAndCrossover[ISLAND_SIZE];
	__shared__ ushort sortedIndices[ISLAND_SIZE];

	FP16 tempChild[3 * 200];

	ushort crossoverIdx;
	FP16 *parent1, *parent2;
	FP16 *localChromosome;
	FP16 a;
	float gauss, mutationProbability;

	curandState localState = state[threadIndex];

	for(int i = 0; i < clusterSize; ++i) {
		population[threadIdx.x * clusterSize + i] = populationGlobal[threadIndex * clusterSize + i];
	}

	int i = 1;

	while(true) {
		localChromosome = &population[threadIdx.x * clusterSize];

		// fitness evaluation
		fitness[threadIdx.x] = calcScoreHeatmapActiveRegion(-1, localChromosome, heatmap_chromosome_boundaries, 
			heatmap_dist, heatmapSize, heatmapDiagonalSize, clusterSize, chromosomeBoundariesSize);

		__syncwarp();

		if(i == ITERATIONS) break;

		// MIGRATION - START
		if(i % 100 == 0) {
			bitonicSort(fitness, sortedIndices, ISLAND_SIZE);
			
			if(threadIdx.x == 0) {
				// write M best chromosomes to neighbouring island on the left
				int islandId = blockIdx.x == 0 ? numberOfIslands - 1 : blockIdx.x - 1;

				for(int j = 0; j < migrationConstant; ++j) {
					for(int k = 0; k < clusterSize; ++k) {
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
			fitness[threadIdx.x] = calcScoreHeatmapActiveRegion(-1, localChromosome, heatmap_chromosome_boundaries, 
				heatmap_dist, heatmapSize, heatmapDiagonalSize, clusterSize, chromosomeBoundariesSize);

			__syncwarp();
		}
		// MIGRATION - END

		// selection and crossover
		crossoverIdx = curand(&localState) % ISLAND_SIZE;
		selectedIndices[threadIdx.x] = __hgt(fitness[threadIdx.x], fitness[crossoverIdx]) ? threadIdx.x : crossoverIdx;
		weightsAndCrossover[threadIdx.x] = __float2half(curand_uniform(&localState));
		__syncwarp();

		crossoverIdx = threadIdx.x % 2 == 0 ? threadIdx.x : threadIdx.x - 1;

		if(__hgt(__float2half(0.7), weightsAndCrossover[crossoverIdx]) {
			a = weightsAndCrossover[crossoverIdx + 1];
			crossoverIdx = threadIdx.x % 2 == 0 ? 1 : -1;  // reuse the variable

			parent1 = &population[selectedIndices[threadIdx.x] * clusterSize];
			parent2 = &population[selectedIndices[threadIdx.x + crossoverIdx] * clusterSize];

			for(int j = 0; j < clusterSize; ++j) {
				// Offspring1 = a * Parent1 + (1- a) * Parent2
				// Offspring2 = (1 - a) * Parent1 + a * Parent2
				tempChild[j] = __hadd(__hmul(a, parent1[j]), __hmul(__hsub(__float2half(1.0), a), parent2[j]));
			}
		}
		__syncwarp();
		
		// mutation and store back to shared memory
		gauss = curand_normal(&localState);
		mutationProbability = curand_uniform(&localState);

		for(int j = 0; j < clusterSize; ++j) {
			localChromosome[j] = tempChild[j];
			if(mutationProbability < 0.05) localChromosome[j] = __hadd(localChromosome[j], __float2half(gauss));
		}
		__syncwarp();
		++i;
	}

	// after we finish, write population back to global memory
	for(int i = 0; i < clusterSize; ++i) {
		populationGlobal[threadIndex * clusterSize + i] = population[threadIdx.x * clusterSize + i];
	}
}

float LooperSolver::ParallelGeneticHeatmap(float step_size) {
	// TODO: initialize gpu data and launch kernel

	const int N = 100;
    const int blocks = 256;  // const int blocks = Settings::numberOfBlocks;
	const int threads = 32;  // const int threads = Settings::numberOfThreads;
    const int numberOfIndividuals = blocks * threads;

	double score_curr;
	double score_density;
	double milestone_score;		// score measured every N steps, used to see whether the solution significantly improves
	int milestone_success = 0;	// calc how many successes there were since last milestone
	string chr; 				// tmp variable to keep track to which chromosome a specific point belongs
    int size = active_region.size();
    int heatmapSize = heatmap_dist.size;

	if (size <= 1) return 0.0;	// there is nothing to do

	// calc initial score
	score_curr = calcScoreHeatmapActiveRegion();	// score from heatmap
	score_density = calcScoreDensity();

	bool * d_hasError;
    bool h_hasError;

    cudaMalloc((void**)&d_hasError, sizeof(bool));
    cudaMemset(d_hasError, 0, sizeof(bool));
    
    struct gpu_settings settings;
    settings.tempJumpScaleHeatmap = Settings::tempJumpScaleHeatmap;
    settings.tempJumpCoefHeatmap = Settings::tempJumpCoefHeatmap;
    settings.use2D = Settings::use2D;
	
	thrust::host_vector<bool> h_clusters_fixed(clusters.size());
    thrust::host_vector<half3> h_clusters_positions(clusters.size()); // just one initial state to copy across
    thrust::host_vector<float> h_scores(numberOfIndividuals);

    for(int i = 0; i < clusters.size(); ++i) {
        h_clusters_fixed[i] = clusters[i].is_fixed;
        h_clusters_positions[i].x = __float2half(clusters[i].pos.x);
        h_clusters_positions[i].y = __float2half(clusters[i].pos.y);
        h_clusters_positions[i].z = __float2half(clusters[i].pos.z);
    }

    thrust::device_vector<int> d_active_region(active_region);
    thrust::device_vector<int> d_heatmap_chromosome_boundaries(heatmap_chromosome_boundaries);
    thrust::device_vector<bool> d_clusters_fixed = h_clusters_fixed;
    thrust::device_vector<half3> d_clusters_positions(numberOfIndividuals * clusters.size());
    thrust::device_vector<float> d_scores(numberOfIndividuals, 0.0);
    thrust::device_vector<float> d_heatmap_dist(heatmapSize * heatmapSize);

    // TODO: make sure it's copying in the correct order
    for(int i = 0; i < heatmapSize; ++i) {
        thrust::copy(heatmap_dist.v[i], heatmap_dist.v[i] + heatmapSize, d_heatmap_dist.begin() + heatmapSize * i);
    }

	for(int i = 0; i < numberOfIndividuals; ++i) {
        thrust::copy(h_clusters_positions.begin(), h_clusters_positions.end(), d_clusters_positions.begin() + i * clusters.size());
    }

	curandState * devStates;
	gpuErrchk(cudaMalloc((void**)&devStates, numberOfIndividuals * sizeof(curandState)));

    // cuRand initialization
	setupKernel<<<blocks, threads>>>(devStates, time(NULL));
	gpuErrchk( cudaDeviceSynchronize() );

	output(2, "initial score: %lf (density=%lf)\n", score_curr, score_density);

	geneticHeatmap<<<blocks, threads>>>(
		devStates, //curandState * state,
		d_clusters_positions, //FP16 * populationGlobal,
		d_heatmap_chromosome_boundaries, //int * heatmap_chromosome_boundaries,
		d_heatmap_dist, //FP16 * heatmap_dist,
		d_scores, //FP16 * scores,
		gpu_settings, //gpu_settings settings,
		blocks, //const int numberOfIslands,
		const int clusterSize, // 3 * active region
		const int migrationConstant,
		const int heatmapSize,
		const int heatmapDiagonalSize,
		const int chromosomeBoundariesSize
	);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	int resultIndex = thrust::min_element(thrust::device, d_scores.begin(), d_scores.end()) - d_scores.begin();

	auto iterStart = d_clusters_positions.begin() + resultIndex * clusters.size();
    auto iterEnd = d_clusters_positions.begin() + resultIndex * clusters.size() + clusters.size();

    for(int i = 0; i < numberOfParallelSimulations; ++i) {
        if(i == resultIndex) continue;
        thrust::copy(iterStart, iterEnd, d_clusters_positions.begin() + i * clusters.size());
	}

	h_scores = d_scores;
	float bestScore = h_scores[resultIndex];
	
	thrust::copy(iterStart, iterEnd, h_clusters_positions.begin());

    for(int i = 0; i < clusters.size(); ++i) {
        clusters[i].pos.x = h_clusters_positions[i].x;
        clusters[i].pos.y = h_clusters_positions[i].y;
        clusters[i].pos.z = h_clusters_positions[i].z;
	}
	
	std::cout << "Score: " << bestScore << std::endl;

	cudaFree(d_hasError);
	cudaFree(devStates);
	return bestScore;
}