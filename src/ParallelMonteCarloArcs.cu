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
__device__ void randomVector(half3& vector, const float& max_size, bool& in2D, curandState* state) {
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

__device__ void subtractFromVector3(half3& destination, half3& value1, half3& value2) {
        destination.x = __hsub(value1.x, value2.x);
        destination.y = __hsub(value1.y, value2.y);
        destination.z = __hsub(value1.z, value2.z);
}

__device__ float magnitude(half3& vector) {
    float x = __half2float(vector.x), y = __half2float(vector.y), z = __half2float(vector.z);
    return sqrtf(x*x + y*y + z*z);
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
    int& warpIdx
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
            if (i == moved || (helper = heatmap_dist[i * heatmapSize + moved]) < 1e-3) continue;	// ignore 0 values
            //no warp divergence as all threads are from the same warp so the warpIdx will evaluate the same
            subtractFromVector3(temp_one, *(clusters_positions + i),
                                          (moved == warpIdx) ? curr_vector : *(clusters_positions + moved));
            helper = magnitude(temp_one) / helper - 1; 
            err += helper * helper;
        }
    }
    return err;
}

__device__ float s_calcScoreHeatmapSingleActiveRegion(
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
    int& warpIdx
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
            if (i == moved || (helper = heatmap_dist[i]) < 1e-3) continue;	// ignore 0 values
            //no warp divergence as all threads are from the same warp so the warpIdx will evaluate the same
            subtractFromVector3(temp_one, *(clusters_positions + i),
                                          (moved == warpIdx) ? curr_vector : *(clusters_positions + moved));
            helper = magnitude(temp_one) / helper - 1; 
            err += helper * helper;
        }
    }
    return err;
}
//this function uses following
//	-> size of the active region
//	-> the whole active region array
//	heatmap_dist
//		-> size
//		-> diagonal size
//		-> entire v 2d array (column at once)
//	heatmap_chromosome_boundaries	
//		-> size
// getChromosomeBoundary
// entire active region array
// all positions from the clusters 
__device__ float calcScoreHeatmapActiveRegion(
    const int moved, 
    half3 * __restrict__ clusters_positions, 
    int * __restrict__ gpu_heatmap_chromosome_boundaries, 
    float * __restrict__ heatmap_dist,
    const int& heatmapSize,
    const int& heatmapDiagonalSize,
    const int& activeRegionSize,
    const int& chromosomeBoundariesSize, 
    half3& curr_vector, 
    int& warpIdx
) {
    float err = 0.0f;
	if (moved == -1) {
        for (int i = 0; i < activeRegionSize; ++i) {
            err += calcScoreHeatmapSingleActiveRegion(i, clusters_positions, gpu_heatmap_chromosome_boundaries, 
                heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize, curr_vector, warpIdx);
        }
        return err;
	}
    return calcScoreHeatmapSingleActiveRegion(moved, clusters_positions, gpu_heatmap_chromosome_boundaries, 
            heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize, curr_vector, warpIdx);
}

// calculate score based on distances between the selected bead and clusters connected with it by arcs
__device__ float calcScoreDistancesActiveRegion(
    const int & cluster_moved, 
    int * __restrict__ active_region, 
    half3 * __restrict__ clusters_positions, 
    float * __restrict__ heatmap_exp_dist_anchor,
    const float & springConstantStretchArcs,
    const float & springConstantSqueezeArcs
) {
	float sc = 0.0f, diff;
	int st = active_region[cluster_moved];
    int n = sizeof(active_region) / sizeof(active_region[0]);
    half3 v;
    
	for (int i = 0; i < n; ++i) {
		if (i == cluster_moved) continue;
        // v = clusters[st].pos - clusters[active_region[i]].pos;
        v.x = clusters_positions[st].x - clusters_positions[active_region[i]].x;
        v.y = clusters_positions[st].y - clusters_positions[active_region[i]].y;
        v.z = clusters_positions[st].z - clusters_positions[active_region[i]].z;

		if (heatmap_exp_dist_anchor[i * n + cluster_moved] < 1e-6) continue;

		diff = (magnitude(v) - heatmap_exp_dist_anchor[i * n + cluster_moved]) / heatmap_exp_dist_anchor[i * n + cluster_moved];
		sc += diff * diff * (diff >= 0.0f ? springConstantStretchArcs : springConstantSqueezeArcs);
	}
	return sc;
}


__global__ void setupKernel(curandState * state, time_t seed) {
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    curand_init(seed, tid, 0, &state[tid]);
}

__global__ void MonteCarloHeatmapKernel(
    curandState * __restrict__ state,
    bool * __restrict__ clusters_fixed,
    half3* __restrict__ clusters_positions,
	int* __restrict__ heatmap_chromosome_boundaries,
	int * __restrict__ milestone_successes,
	float * __restrict__ scores,
	float T,
    const float score,
	gpu_settings settings,
	float * __restrict__ heatmap_dist,
    const float step_size,
    const int numberOfClusters,
    const int activeRegionSize,
    const int heatmapSize,
    const int heatmapDiagonalSize,
    const int chromosomeBoundariesSize,
    bool * __restrict__ error
) {
    int iterations = 0;
    //TODO bring these back or delete them
    //float milestone_score = score;
    //int laneIdx = threadIdx.x & (warpSize - 1); // we know that warp size is 32 so we can change modulo for bitwise and
    // this could be given as argument (e.g. 2 warps per el then warpIdx = threadIndex/(warpSize * 2))
    // should we use a register? operation may be costly but it's a few cycles costly while memory operation is 
    // few hundred cycles costly
    extern __shared__ float s_heatmap_dist[];
    
    int threadIndex = blockDim.x * blockIdx.x + threadIdx.x;
    int warpIdx = blockIdx.x % activeRegionSize;
    //int warpIdx = (threadIndex / warpSize) % activeRegionSize;
    half3 curr_vector;
    half3 displacement;
	float score_curr = score;
    float score_prev = score_curr;
    __shared__ unsigned int mutex;
    #if __CUDA_ARCH >= 800
        int i_winner;
        int i_score;
    #else
        float f_winner;
    #endif
    curandState localState = state[threadIndex];
    for(int i = 0; i <= activeRegionSize / blockDim.x; ++i) {
        if(threadIdx.x + i * blockDim.x < activeRegionSize)
            s_heatmap_dist[threadIdx.x + i * blockDim.x] = heatmap_dist[activeRegionSize * (threadIdx.x + i * blockDim.x) + warpIdx];
    }
    __syncthreads();
    // MOVE ALL RELEVANT POSITIONS INTO THE LOCAL MEM
    while(true) {
        if(threadIndex == 0) {
                  printf("GPU CALC FUNCTION VERIFICATION\n CPU SCORE: %f GPU SCORE %f \n", score, calcScoreHeatmapActiveRegion(-1, clusters_positions, heatmap_chromosome_boundaries, 
                       heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize, *(clusters_positions + warpIdx), warpIdx));
                }
        curr_vector = clusters_positions[warpIdx];
        #define N 512
        #pragma unroll
        for(int i = 0; i < N; ++i) {
            if (clusters_fixed[warpIdx]) *error = true;
            if(*error) return;
            randomVector(displacement, step_size , settings.use2D, &localState);
            addToVector(curr_vector, displacement);
            // TODO could there be a way to recompute only the change instead of the whole score? pen and paper
            score_curr = calcScoreHeatmapSingleActiveRegion(warpIdx, clusters_positions, heatmap_chromosome_boundaries, 
                    heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize, curr_vector, warpIdx);
            if((score_curr <= score_prev)
                // TODO can this be precomputed in some part (exp is costly)?
               || ( T > 0.0f && withChance(settings.tempJumpScaleHeatmap * expf(-settings.tempJumpCoefHeatmap * (score_curr * (1/score_prev)) * (1/T)), &localState))) {
                // not thread safe, does it matter?
                ++milestone_successes[0];
                score_prev = score_curr;
                continue;
            }
            score_curr = score_prev;
            subtractFromVector(curr_vector, displacement);
        }
        T *= 0.9999;
        iterations += N;
        //find the best move
        __syncwarp();
        #define FULL_MASK 0xffffffff
        #if __CUDA_ARCH__ >= 800
            i_score = (int)score_curr; 
            i_winner = __reduce_min_sync(mask, score_red);
            if(i_winner == i_score) {
        #else
            f_winner = score_curr;
            #pragma unroll 5
            for (int offset = 16; offset > 0; offset /= 2)
                f_winner = fminf(f_winner, __shfl_down_sync(FULL_MASK, f_winner, offset));
            //propagate and check who's the lucky one
            __syncwarp();
            __shfl_sync(FULL_MASK, f_winner, 0);
            if(f_winner == score_curr) {
        #endif
            // TODO doubling of the same work here -> memoization or specialized function for this point
            //mutex_lock(&mutex);
            if(calcScoreHeatmapSingleActiveRegion(warpIdx, clusters_positions, heatmap_chromosome_boundaries, 
                heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize, curr_vector, warpIdx) <
                calcScoreHeatmapSingleActiveRegion(warpIdx, clusters_positions, heatmap_chromosome_boundaries, 
                heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize, *(clusters_positions + warpIdx), warpIdx)) {
                    clusters_positions[warpIdx] = curr_vector;
                }
            //mutex_unlock(&mutex);
        }
        // TODO THREADFENCE OR THREADFENCE BLOCK
        __threadfence();
        #if __CUDA_ARCH__ >= 800
            score_curr = i_winner;
        #else
            score_curr = f_winner;
        #endif
        score_prev = score_curr;
        // check if we should stop
        // make sure that Settings::MCstopConditionSteps is divisible by N
        // 32767 = 32768 - 1, which is a power of 2
        if (iterations % 16384 == 0) {
            //if ((score_curr > Settings::MCstopConditionImprovementHeatmap * milestone_score &&
             //       milestone_success < Settings::MCstopConditionMinSuccessesHeatmap) || score_curr < 1e-6) {
                
                scores[0] = score_curr;
                break;
            //}
            //else {
            //    milestone_score = score_curr;
           // }
            
        }
    }
    state[threadIndex] = localState;
}




float LooperSolver::ParallelMonteCarloHeatmap(float step_size) {
    const int blocks = active_region.size();
    const int threads = 256;
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
    thrust::host_vector<float> h_clusters_positions(active_region.size() * 3); // just one initial state to copy across
    thrust::host_vector<float> h_scores(1);
    thrust::host_vector<int> h_milestone_successes(1);

    for(int i = 0; i < active_region.size(); ++i) {
        h_clusters_fixed[i] = clusters[i].is_fixed;
        #pragma unroll
        for(int j =0; j < 3; ++j) {
            h_clusters_positions[i*3 + j] = clusters[active_region[i]].pos[j];
        }
    }
    
    thrust::device_vector<int> d_heatmap_chromosome_boundaries(heatmap_chromosome_boundaries);
    thrust::device_vector<bool> d_clusters_fixed = h_clusters_fixed;
    thrust::host_vector<half3> h_clusters_positions_half(active_region.size());
    thrust::device_vector<half3> d_clusters_positions(active_region.size());
	thrust::device_vector<int> d_milestone_successes(1, 0); 
    thrust::device_vector<float> d_scores(1, 0.0f);
    thrust::device_vector<float> d_heatmap_dist(heatmapSize * heatmapSize);
    //TODO JUST FOR TEST PURPOSES

    // TODO: make sure it's copying in the correct order
    for(int i = 0; i < heatmapSize; ++i) {
        thrust::copy(heatmap_dist.v[i], heatmap_dist.v[i] + heatmapSize, d_heatmap_dist.begin() + heatmapSize * i);
    }
    for(int i = 0; i < active_region.size(); ++i) {
        h_clusters_positions_half[i].x = __float2half(clusters[active_region[i]].pos.x);
        h_clusters_positions_half[i].y = __float2half(clusters[active_region[i]].pos.y);
        h_clusters_positions_half[i].z = __float2half(clusters[active_region[i]].pos.z);
    }
    //MOVE NON CHANGEABLE THINGS TO CONSTANT MEMORY
    //TODO depending on the result of printing we can get away with that or not
    d_clusters_positions = h_clusters_positions_half;
	curandState * devStates;
	gpuErrchk(cudaMalloc((void**)&devStates, threads * blocks * sizeof(curandState)));

    // cuRand initialization
	setupKernel<<< blocks, threads>>>(devStates, time(NULL));
	gpuErrchk( cudaDeviceSynchronize() );

	output(2, "initial score: %lf (density=%lf)\n", score_curr, score_density);
    std::cout << "CHROMOSOME BOUNDARIES SIZE IS : " << heatmap_chromosome_boundaries.size() << std::endl;
    MonteCarloHeatmapKernel<<<blocks, threads, active_region.size() * sizeof(float)>>>(
        devStates,
        thrust::raw_pointer_cast(d_clusters_fixed.data()),
        thrust::raw_pointer_cast(d_clusters_positions.data()),
        thrust::raw_pointer_cast(d_heatmap_chromosome_boundaries.data()),
        thrust::raw_pointer_cast(d_milestone_successes.data()),
        thrust::raw_pointer_cast(d_scores.data()),
        T,
        score_curr,
        settings,
        thrust::raw_pointer_cast(d_heatmap_dist.data()),
        // TODO BIG STEP SIZE FOR BIG CHROMOSOMES SMALL FOR SMALL
        0.1f * step_size,
        clusters.size(),
        active_region.size(),
        heatmapSize,
        heatmap_dist.diagonal_size,
        heatmap_chromosome_boundaries.size(),
        d_hasError
    );

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
    h_clusters_positions_half = d_clusters_positions;
    #pragma unroll 58
    for(int i = 0; i < active_region.size(); ++i) {
        clusters[active_region[i]].pos.x = __half2float(h_clusters_positions_half[i].x);
        clusters[active_region[i]].pos.y = __half2float(h_clusters_positions_half[i].y);
        clusters[active_region[i]].pos.z = __half2float(h_clusters_positions_half[i].z);
    }
    printf("===========================================================");
    printf("CPU SCORE IS %f \n", calcScoreHeatmapActiveRegion());
	// TODO: do we need more cudaFree?
	cudaFree(d_hasError);
    cudaFree(devStates);
	return score_curr;
}
