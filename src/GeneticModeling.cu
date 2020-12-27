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
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ void bitonicSort(float *dev_values, ushort *indices, int islandSize) {
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
					if (dev_values[i] > dev_values[ixj]) {
						/* exchange(i,ixj); */
						float temp = dev_values[i];
						dev_values[i] = dev_values[ixj];
						dev_values[ixj] = temp;
						indices[i] = ixj;
						indices[ixj] = i;
					}
				}
				if ((i&k) != 0) {
					/* Sort descending */
					if (dev_values[i] < dev_values[ixj]) {
						/* exchange(i,ixj); */
						float temp = dev_values[i];
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

// ACTION ITEMS:
// 1. use FP16 to allow for 2x larger population in a block
// 2. try doing crossover on one thread from each pair only (less barriers, but one thread idle)
// 		is it even possible?
// 3. optimize migration by utilizing more than one thread for data transfer
// 4. Since we probably won't be able to fit more than 32 elements per block, we can use warp-level operations
// 5. Can we do better with the calcScoreHeatmapActiveRegion function?
// 6. We could utilize Tensor Cores for offspring calculations
__global__ void geneticHeatmap(
	curandState * state,
	float * populationGlobal,
	int * heatmap_chromosome_boundaries,
	int * milestone_successes,
	float * heatmap_dist,
	float * scores,
	gpu_settings settings,
    const int numberOfIslands,
	const int clusterSize, // 3 * active region
	const int migrationConstant,
    const int heatmapSize,
    const int heatmapDiagonalSize,
    const int chromosomeBoundariesSize
) {
	int threadIndex = blockDim.x * blockIdx.x + threadIdx.x;
	if(threadIdx.x >= ISLAND_SIZE) return;

	extern __shared__ float population[];
	__shared__ float fitness[ISLAND_SIZE];
	__shared__ ushort selectedIndices[ISLAND_SIZE];
	__shared__ float weightsAndCrossover[ISLAND_SIZE];
	__shared__ ushort sortedIndices[ISLAND_SIZE];

	float tempChild[3 * 200];

	ushort crossoverIdx;
	float *parent1, *parent2;
	float *localChromosome;
	float a, gauss, mutationProbability;

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
		selectedIndices[threadIdx.x] = fitness[threadIdx.x] > fitness[crossoverIdx] ? threadIdx.x : crossoverIdx;
		weightsAndCrossover[threadIdx.x] = curand_uniform(&localState);
		__syncwarp();

		crossoverIdx = threadIdx.x % 2 == 0 ? threadIdx.x : threadIdx.x - 1;

		if(weightsAndCrossover[crossoverIdx] < 0.7) {
			a = weightsAndCrossover[crossoverIdx + 1];
			crossoverIdx = threadIdx.x % 2 == 0 ? 1 : -1;  // reuse the variable

			parent1 = &population[selectedIndices[threadIdx.x] * clusterSize];
			parent2 = &population[selectedIndices[threadIdx.x + crossoverIdx] * clusterSize];

			for(int j = 0; j < clusterSize; ++j) {
				// Offspring1 = a * Parent1 + (1- a) * Parent2
				// Offspring2 = (1 - a) * Parent1 + a * Parent2
				tempChild[j] = a * parent1[j] + (1.0 - a) * parent2[j];
			}
		}
		__syncwarp();
		
		// mutation and store back to shared memory
		gauss = curand_normal(&localState);
		mutationProbability = curand_uniform(&localState);

		for(int j = 0; j < clusterSize; ++j) {
			localChromosome[j] = tempChild[j];
			if(mutationProbability < 0.05) localChromosome[j] += gauss;
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
	
	
}