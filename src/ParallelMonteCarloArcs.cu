#include <thrust/fill.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <cuda.h>
#include <cuda_fp16.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <float.h>
#include <time.h>

#include "../include/LooperSolver.h"


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

struct gpu_settings {
	float densityInfluence;
	float densityScale;
	float densityWeight;
	float maxTempHeatmapDensity;
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
};

struct __align__(8) half3{
    __half x;
    __half y;
    __half z;
};

__device__ void mutex_lock(unsigned int *mutex) {
    unsigned int ns = 10;
    while (atomicCAS(mutex, 0, 1) == 1) {
        __nanosleep(ns);
        if (ns < 256) {
        ns *= 2;
        }
    }
}

__device__ void mutex_unlock(unsigned int *mutex) {
    atomicExch(mutex, 0);
}

__device__ int random(int range, curandState * state) {
	return curand(state) % range;
}

//todo const references
__device__ __half random(const float& range, bool negative, curandState * state) {
	if (negative) return __float2half((2.0f * curand_uniform(state) - 1.0f) * range);
	return __float2half(range * curand_uniform(state));
}


//todo const references
__device__ bool withChance(const float& chance, curandState* state) {
	return curand_uniform(state) < chance;
}

//todo const references
__device__ void randomVector(half3& vector, const float& max_size, const bool& in2D, curandState* state) {
    vector.x = random(max_size, true, state);
    vector.y = random(max_size, true, state);
    vector.z = in2D ? __float2half(0.0f) : random(max_size, true, state);
}

__device__ void addToVector(half3& destination, half3& value) {
    destination.x = __hadd(destination.x, value.x);
    destination.y = __hadd(destination.y, value.y);
    destination.z = __hadd(destination.z, value.z);
}

__device__ void subtractFromVector(half3& destination, half3& value) {
        destination.x = __hsub(destination.x, value.x);
        destination.y = __hsub(destination.y, value.y);
        destination.z = __hsub(destination.z, value.z);
}

__device__ void add3Vectors(half3& dest, const half3& a, const half3& b, const half3& c) {
        dest.x = __hadd(c.x, __hadd(a.x, b.x));
        dest.y = __hadd(c.y, __hadd(a.y, b.y));
        dest.z = __hadd(c.z, __hadd(a.z, b.z));
}

__device__ void subtractFromVector3(half3& destination, half3& value1, half3& value2) {
        destination.x = __hsub(value1.x, value2.x);
        destination.y = __hsub(value1.y, value2.y);
        destination.z = __hsub(value1.z, value2.z);
}

__device__ half3 subtractAndReturnVector(const half3& a, const half3& b) {
        half3 ret;
        ret.x = __hsub(a.x, b.x);
        ret.y = __hsub(a.y, b.y);
        ret.z = __hsub(a.z, b.z);
        return ret;
}
__device__ half3 multiplyVector(const half3& v, const float& fv) {
        half3 ret;
        ret.x = __hmul(v.x, __float2half(fv));
        ret.y = __hmul(v.y, __float2half(fv));
        ret.z = __hmul(v.z, __float2half(fv));
        return ret;
}

__device__ float magnitude(half3& vector) {
    float x = __half2float(vector.x), y = __half2float(vector.y), z = __half2float(vector.z);
    return sqrtf(x*x + y * y + z * z);

}

// TODO make it inline
__device__
float warpReduceMin(float val) {
    __syncwarp();
    unsigned mask = __activemask(); 
    #if __CUDA_ARCH__ >= 800
        val = __reduce_min_sync(mask, val);
    #else
        #pragma unroll 5
        for (int offset = 16; offset > 0; offset /= 2)
            val = fminf(val, __shfl_down_sync(mask, val, offset));
        //propagate and check who's the lucky one
        __syncwarp();
        val = __shfl_sync(mask, val, 0);
    #endif
    return val;
}

// FOR AT MOST 32 WARPS PER BLOCK
// TODO inline it
__device__
float blockReduceMin(float val) {
  static __shared__ float shared[32]; // Shared mem for 32 partial sums
  int lane = threadIdx.x % warpSize;
  int wid = threadIdx.x / warpSize;

  val = warpReduceMin(val);     // Each warp performs partial reduction

  if (lane==0) shared[wid]=val; 

  __syncthreads(); 
  val = (threadIdx.x <= blockDim.x / warpSize) ? shared[lane] : FLT_MAX;

  if (wid==0) { shared[lane] = warpReduceMin(val); }//Final reduce within first warp
  __syncthreads();
  return shared[0];
}

__inline__ __device__
float warpReduceAdd(float val) {
    __syncwarp();
    unsigned mask = __activemask(); 
    #if __CUDA_ARCH__ >= 800
        val = __reduce_add_sync(mask, val);
    #else
        #pragma unroll 5
        for (int offset = 16; offset > 0; offset /= 2)
            val += __shfl_down_sync(mask, val, offset);
        __syncwarp();
        val = __shfl_sync(mask, val, 0);
    #endif
    return val;
}

// FOR AT MOST 32 WARPS PER BLOCK
__inline__ __device__
float blockReduceAdd(float val) {
  static __shared__ float shared[32];
  int lane = threadIdx.x % warpSize;
  int wid = threadIdx.x / warpSize;
  val = warpReduceAdd(val);   
  if (lane==0) shared[wid]=val;
  __syncthreads();             
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid==0) { val = warpReduceAdd(val); }//Final reduce within first warp
  __syncthreads();
  return val;
}

//this function uses following external funcitonality / data
// -> size of heatmap_chromosome_boundaries
// -> the entire heatmap_chromosome_boundaries
__device__ void getChromosomeHeatmapBoundary(const int p, int &start, int &end, int* gpu_heatmap_chromosome_boundaries, const int& chromosomeBoundariesSize) {
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
    const int moved, 
    half3* __restrict__ clusters_positions, 
    int* __restrict__ gpu_heatmap_chromosome_boundaries, 
    float * __restrict__ heatmap_dist,
    const int& heatmapSize,
    const int& heatmapDiagonalSize,
    const int& activeRegionSize,
    const int& chromosomeBoundariesSize, 
    //TODO this reference may be a bad idea 
    half3& curr_vector,
    const int& warpIdx
) {
    double err = 0.0;
    float helper;
    int st = 0, end = activeRegionSize - 1;
    if (heatmapSize != activeRegionSize) { printf("heatmap sizes mismatch, dist size=%d, active region=%d", heatmapSize, activeRegionSize); return 0.0; }
    half3 temp_one;
    
    // TODO can we precompute?
    // array with start and end 
    if(chromosomeBoundariesSize > 0) {
        getChromosomeHeatmapBoundary(moved, st, end, gpu_heatmap_chromosome_boundaries, chromosomeBoundariesSize);
    }

    for (int i = st; i <= end; ++i) {
        if (abs(i-moved) >= heatmapDiagonalSize) {
            if (i == moved || (helper = heatmap_dist[i * heatmapSize + moved]) < 1e-6) continue;	// ignore 0 values
            //no warp divergence as all threads are from the same warp so the warpIdx will evaluate the same
            subtractFromVector3(temp_one, *(clusters_positions + i),
                                          (moved == warpIdx) ? curr_vector : *(clusters_positions + moved));
            helper = magnitude(temp_one) * (1 / helper) - 1;
            err += helper * helper;
        }
    }
    return err;
}

__global__ void setupKernel(curandState * state, time_t seed) {
    int posid = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
    curand_init(seed, posid, 0, &state[posid]);
}


__global__ void initializationKernel(
    curandState * __restrict__ state,
    half3 * positions,
    half3 * velocities
) {
    int posID = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
    // update the initially given positon 
    half3 rnd;
    randomVector(rnd, 30.0f, false, (state + posID));
    addToVector(positions[posID], rnd);
    randomVector(velocities[posID], 60.0f, false, (state + posID));
}

__global__ void positionUpdateKernel(
    curandState * __restrict__ state,
    half3* positions, 
    half3* local_best_positions,
    half3* velocities,
    float* fitnesses,
    float* local_best_fitnesses,
    float* global_best_fitnesses, 
    half3* global_best_positions,
    int generation,
    int maxGenerations
) {
   // number of thread blocks equal to the number of particles
   // each thread updates the position of one dim of one particle - block per part
   // global best, local best need to be passed as args 
    #define WMAX 0.9
    #define WMIN 0.4
    #define C1 2.05
    #define C2 2.05
    int posID = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
    float inertia = WMAX - ((WMAX - WMIN) / maxGenerations) * generation;
    half3 temp1 = multiplyVector(velocities[posID], inertia);
    half3 temp2 = multiplyVector(subtractAndReturnVector(global_best_positions[blockIdx.x], positions[posID]), curand_uniform(state + posID) * C1);
    half3 temp3 = multiplyVector(subtractAndReturnVector(local_best_positions[posID], positions[posID]), curand_uniform(state + posID) * C2);
    add3Vectors(velocities[posID], temp1, temp2, temp3);
    addToVector(positions[posID], velocities[posID]);
}

__global__ void fitnessKernel(
    half3* __restrict__ positions,
    float* fitnesses,
    bool * __restrict__ clusters_fixed,
	int* __restrict__ heatmap_chromosome_boundaries,
	gpu_settings settings,
	float * __restrict__ heatmap_dist,
    const int activeRegionSize,
    const int heatmapSize,
    const int heatmapDiagonalSize,
    const int chromosomeBoundariesSize
) {
    // whole block is responsible for one particle 
    // set of blocks in a swarm is responsible for a swarm
    int posID = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
    float local_score = 0;
    if(threadIdx.x < activeRegionSize)
        local_score = calcScoreHeatmapSingleActiveRegion(threadIdx.x, positions + posID - threadIdx.x, heatmap_chromosome_boundaries, 
        heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize, positions[posID], threadIdx.x);
    __syncthreads();
    local_score = blockReduceAdd(local_score);
    if(threadIdx.x == 0) {
        fitnesses[blockIdx.x * gridDim.y + blockIdx.y] = local_score;
    }
}

// TODO This one is easily more parallelizable
__global__ void bestUpdateKernel(
    float* fitnesses,
    float* local_best_fitnesses,
    float* global_best_fitnesses,
    half3* positions,
    half3* local_best_positions,
    half3* global_best_positions,
    int numDimensions
){
    // one thread per particle
    int particleID = blockIdx.x * blockDim.x + threadIdx.x;
    // for each particle check if current fitness < local_best_fitness
    if(local_best_fitnesses[particleID] >= fitnesses[particleID]) {
        local_best_fitnesses[particleID] = fitnesses[particleID];
        // easily parallelizable
        for(int i = 0; i < numDimensions; ++i) {
            local_best_positions[particleID * numDimensions + i] = positions[(blockIdx.x * blockDim.x + threadIdx.x) * numDimensions + i];
        }
    }

    __syncthreads();
    // reduce local_best_positions into global best for swarm
    float local_best = local_best_fitnesses[particleID];
    __syncwarp();
    float total_local_best = blockReduceMin(local_best);
    __syncthreads();
    if(total_local_best < global_best_fitnesses[blockIdx.x] && local_best == total_local_best) {
        global_best_fitnesses[blockIdx.x] = local_best;
        for(int i = 0; i < numDimensions; ++i) {
            global_best_positions[blockIdx.x * numDimensions + i] = local_best_positions[particleID * numDimensions + i];
        }
    }
}

float LooperSolver::ParallelMonteCarloHeatmap(float step_size) {
    const auto numSwarms = 32;
    const auto numParticles = 320;
    const auto numDimensions = active_region.size();
    // const int blocks = Settings::numberOfBlocks;
    // const int threads = Settings::numberOfThreads;

	double T = Settings::maxTempHeatmap;		// set current temperature to max

	double score_curr;
	double score_density;
	int milestone_success = 0;	// calc how many successes there were since last milestone
	string chr; 				// tmp variable to keep track to which chromosome a specific point belongs
    int size = active_region.size();
    int heatmapSize = heatmap_dist.size;

	if (size <= 1) return 0.0f;	// there is nothing to do

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
	
	thrust::host_vector<bool> h_clusters_fixed(active_region.size());
    thrust::host_vector<float> h_scores(1);
    thrust::host_vector<int> h_milestone_successes(1);
    thrust::host_vector<float> h_fitnesses(numSwarms * numParticles);
    thrust::host_vector<float> h_best_fitnesses(numParticles * numSwarms);
    thrust::host_vector<float> h_global_best_fitnesses(numSwarms);
    thrust::host_vector<half3> h_best_positions(numParticles * numDimensions);
    thrust::host_vector<half3> h_global_best_positions(numParticles * numDimensions);
    thrust::host_vector<half3> h_clusters_positions(numSwarms * numParticles * numDimensions);
    thrust::host_vector<half3> h_velocities(numSwarms * numParticles * numDimensions);
    
    thrust::device_vector<float> d_fitnesses(numSwarms * numParticles);
    thrust::fill(d_fitnesses.begin(), d_fitnesses.end(), FLT_MAX);
    thrust::device_vector<float> d_best_fitnesses(numParticles * numSwarms);
    thrust::fill(d_best_fitnesses.begin(), d_best_fitnesses.end(), FLT_MAX);
    thrust::device_vector<float> d_global_best_fitnesses(numSwarms);
    thrust::fill(d_global_best_fitnesses.begin(), d_global_best_fitnesses.end(), FLT_MAX);
    thrust::device_vector<half3> d_global_best_positions(numParticles * numDimensions);
    thrust::device_vector<half3> d_best_positions(numParticles * numDimensions);
    thrust::device_vector<half3> d_clusters_positions(numSwarms * numParticles * numDimensions);
    thrust::device_vector<half3> d_velocities(numSwarms * numParticles * numDimensions);

    for(int j = 0; j < numSwarms * numParticles; ++j){
        for(int i = 0; i < active_region.size(); ++i) {
            h_clusters_positions[j * active_region.size() + i].x = clusters[active_region[i]].pos.x;
            h_clusters_positions[j * active_region.size() + i].y = clusters[active_region[i]].pos.y;
            h_clusters_positions[j * active_region.size() + i].z = clusters[active_region[i]].pos.z;
        }
    }
    d_clusters_positions = h_clusters_positions;
    thrust::device_vector<int> d_heatmap_chromosome_boundaries(heatmap_chromosome_boundaries);
    thrust::device_vector<bool> d_clusters_fixed = h_clusters_fixed;
	thrust::device_vector<int> d_milestone_successes(1, 0); 
    thrust::device_vector<float> d_scores(1, 0.0f);
    thrust::device_vector<float> d_heatmap_dist(heatmapSize * heatmapSize);
    for(int i = 0; i < heatmapSize; ++i) {
        thrust::copy(heatmap_dist.v[i], heatmap_dist.v[i] + heatmapSize, d_heatmap_dist.begin() + heatmapSize * i);
    }
    //MOVE NON CHANGEABLE THINGS TO CONSTANT MEMORY
	curandState * devStates;
	gpuErrchk(cudaMalloc((void**)&devStates, numSwarms * numParticles * numDimensions * sizeof(curandState)));

    // cuRand initialization
	setupKernel<<< dim3(numSwarms, numParticles), numDimensions>>>(devStates, time(NULL));
	gpuErrchk( cudaDeviceSynchronize() );

    output(2, "initial score: %lf (density=%lf)\n", score_curr, score_density);
    initializationKernel<<<dim3(numSwarms, numParticles), numDimensions>>>(
        devStates,
        thrust::raw_pointer_cast(d_clusters_positions.data()),
        thrust::raw_pointer_cast(d_velocities.data())
    );
    fitnessKernel<<<dim3(numSwarms, numParticles), numDimensions>>>(
            thrust::raw_pointer_cast(d_clusters_positions.data()),
            thrust::raw_pointer_cast(d_fitnesses.data()),
            thrust::raw_pointer_cast(d_clusters_fixed.data()),
            thrust::raw_pointer_cast(d_heatmap_chromosome_boundaries.data()),
            settings,
            thrust::raw_pointer_cast(d_heatmap_dist.data()),
            active_region.size(),
            heatmapSize,
            heatmap_dist.diagonal_size,
            heatmap_chromosome_boundaries.size()
    );
    bestUpdateKernel<<<numSwarms, numParticles>>>(
            thrust::raw_pointer_cast(d_fitnesses.data()),
            thrust::raw_pointer_cast(d_best_fitnesses.data()),
            thrust::raw_pointer_cast(d_global_best_fitnesses.data()),
            thrust::raw_pointer_cast(d_clusters_positions.data()),
            thrust::raw_pointer_cast(d_best_positions.data()),
            thrust::raw_pointer_cast(d_global_best_positions.data()),
            numDimensions
    );
    #define NUM_GENERATIONS 20000
    for(auto generation = 0; generation < NUM_GENERATIONS; ++generation) {
        gpuErrchk( cudaDeviceSynchronize() );
        positionUpdateKernel<<<dim3(numSwarms, numParticles), numDimensions>>>(
                devStates,
                thrust::raw_pointer_cast(d_clusters_positions.data()),
                thrust::raw_pointer_cast(d_best_positions.data()),
                thrust::raw_pointer_cast(d_velocities.data()),
                thrust::raw_pointer_cast(d_fitnesses.data()),
                thrust::raw_pointer_cast(d_best_fitnesses.data()),
                thrust::raw_pointer_cast(d_global_best_fitnesses.data()),
                thrust::raw_pointer_cast(d_global_best_positions.data()),
                generation,
                NUM_GENERATIONS
        );
        fitnessKernel<<<dim3(numSwarms, numParticles), numDimensions + 32 - numDimensions % 32>>>(
                thrust::raw_pointer_cast(d_clusters_positions.data()),
                thrust::raw_pointer_cast(d_fitnesses.data()),
                thrust::raw_pointer_cast(d_clusters_fixed.data()),
                thrust::raw_pointer_cast(d_heatmap_chromosome_boundaries.data()),
                settings,
                thrust::raw_pointer_cast(d_heatmap_dist.data()),
                active_region.size(),
                heatmapSize,
                heatmap_dist.diagonal_size,
                heatmap_chromosome_boundaries.size()
        );
        bestUpdateKernel<<<numSwarms, numParticles>>>(
                thrust::raw_pointer_cast(d_fitnesses.data()),
                thrust::raw_pointer_cast(d_best_fitnesses.data()),
                thrust::raw_pointer_cast(d_global_best_fitnesses.data()),
                thrust::raw_pointer_cast(d_clusters_positions.data()),
                thrust::raw_pointer_cast(d_best_positions.data()),
                thrust::raw_pointer_cast(d_global_best_positions.data()),
                numDimensions
        );
    }
    gpuErrchk(cudaDeviceSynchronize());
    h_best_fitnesses = d_best_fitnesses;
    for(int i = 0; i < numSwarms; ++i) {
        printf(" %f ", h_best_fitnesses[i]);
    }
    printf("\n");
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
    cudaMemcpy(&h_hasError, d_hasError, sizeof(bool), cudaMemcpyDeviceToHost);
    if(h_hasError) error("cluster fixed during arcs!\n");
    //h_scores = d_scores;
    h_milestone_successes = d_milestone_successes;
    score_curr = h_scores[0];
    printf("FINAL SCORE IS %f", score_curr);
    milestone_success = h_milestone_successes[0];
    // TODO may require a copying kernel
    //h_clusters_positions = d_clusters_positions;
    h_clusters_positions = d_clusters_positions;
    #pragma unroll 58
    for(int i = 0; i < active_region.size(); ++i) {
        clusters[active_region[i]].pos.x = __half2float(h_clusters_positions[i].x);
        clusters[active_region[i]].pos.y = __half2float(h_clusters_positions[i].y);
        clusters[active_region[i]].pos.z = __half2float(h_clusters_positions[i].z);
    }
    printf("===========================================================");
    printf("CPU SCORE IS %f \n", calcScoreHeatmapActiveRegion());
	// TODO: do we need more cudaFree?
	cudaFree(d_hasError);
    cudaFree(devStates);
	return score_curr;
}
