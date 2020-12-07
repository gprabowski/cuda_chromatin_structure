#include <thrust/fill.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>

#include <cuda.h>
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
};

__device__ int random(int range, curandState * state) {
	return curand(state) % range;
}

__device__ float random(float range, bool negative, curandState * state) {
	if (negative) return (2.0f * curand_uniform(state) - 1.0f) * range;
	return range * curand_uniform(state);
}

__device__ bool withChance(float chance, curandState * state) {
	return curand_uniform(state) < chance;
}

__device__ void randomVector(float *vector, float max_size, bool in2D, curandState * state) {
	if (in2D) {
        vector[0] = random(max_size, true, state);
        vector[1] = random(max_size, true, state);
        vector[2] = 0.0;
    } else {
        for(int i = 0; i < 3; ++i) {
            vector[i] = random(max_size, true, state);
        }
    }
}

__device__ void addToVector(float *destination, float *value) {
    #pragma unroll
    for(int i = 0; i < 3; ++i) {
        destination[i] += value[i];
    }
}

__device__ void subtractFromVector(float *destination, float *value) {
    for(int i = 0; i < 3; ++i) {
        *destination -= value[i];
        ++destination;
    }
}

__device__ void subtractFromVector3(float *destination, const float *value1, const float *value2) {
	#pragma unroll
    for(int i = 0; i < 3; ++i) {
        destination[i] = value1[i] - value2[i];
    }
}

__device__ float magnitude(float * vector) {
    int sum = 0;

    for(int i = 0; i < 3; ++i) {
        sum += vector[i] * vector[i];
    }
    return sqrtf(sum);
}

//this function uses following external funcitonality / data
// -> size of heatmap_chromosome_boundaries
// -> the entire heatmap_chromosome_boundaries
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

__device__ double calcScoreHeatmapSingleActiveRegion(
    int moved, 
    float* clusters_positions, 
    int* gpu_heatmap_chromosome_boundaries, 
    float * heatmap_dist,
    int heatmapSize,
    int heatmapDiagonalSize,
    int activeRegionSize,
    int chromosomeBoundariesSize
) {
    double err = 0.0, cerr;
	float temp[3];
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
            const float* temp_one = clusters_positions + i * 3;
            const float* temp_two = clusters_positions + moved * 3;
            subtractFromVector3(temp, temp_one, temp_two);
            d = magnitude(temp);
            cerr = (d - heatmap_dist[i * heatmapSize + moved]) / heatmap_dist[i * heatmapSize + moved];
            err += cerr * cerr;
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
__device__ double calcScoreHeatmapActiveRegion(
    int moved, 
    float* clusters_positions, 
    int* gpu_heatmap_chromosome_boundaries, 
    float * heatmap_dist,
    int heatmapSize,
    int heatmapDiagonalSize,
    int activeRegionSize,
    int chromosomeBoundariesSize
) {
    double err = 0.0;
    size_t n = activeRegionSize;
	if (moved == -1) {
        for (size_t i = 0; i < n; ++i) {
            err += calcScoreHeatmapSingleActiveRegion(i, clusters_positions, gpu_heatmap_chromosome_boundaries, 
                heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize);
        }
	}
	else {
        err = calcScoreHeatmapSingleActiveRegion(moved, clusters_positions, gpu_heatmap_chromosome_boundaries, 
                heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize);
	}
	return err;
}

// calculate score based on distances between the selected bead and clusters connected with it by arcs
__device__ float calcScoreDistancesActiveRegion(
    const int & cluster_moved, 
    int * active_region, 
    float * clusters_positions, 
    float * heatmap_exp_dist_anchor,
    const float & springConstantStretchArcs,
    const float & springConstantSqueezeArcs
) {
	float sc = 0.0, diff;
	int st = active_region[cluster_moved];
    int n = sizeof(active_region) / sizeof(active_region[0]);
    float v[3];
    
	for (int i = 0; i < n; ++i) {
		if (i == cluster_moved) continue;

        // v = clusters[st].pos - clusters[active_region[i]].pos;
        for(int k = 0; k < 3; ++k) {
            v[k] = clusters_positions[st * 3 + k] - clusters_positions[active_region[i] * 3 + k];
        }

        if (heatmap_exp_dist_anchor[i * n + cluster_moved] < 0.0f) {
			//	sc += 1.0f / v.length();
			//sc += 1.0f;
			continue;
		}

		if (heatmap_exp_dist_anchor[i * n + cluster_moved] < 1e-6) continue;	// anchors not connected by an arc are denoted by -1

		diff = (magnitude(v) - heatmap_exp_dist_anchor[i * n + cluster_moved]) / heatmap_exp_dist_anchor[i * n + cluster_moved];
		//printf("%d %d (%d %d)  -> (%lf %lf) %lf\n", i, cluster_moved, st, active_region[i], v.length(), heatmap_exp_dist_anchor.v[i][cluster_moved], diff);
		sc += diff * diff * (diff >= 0.0f ? springConstantStretchArcs : springConstantSqueezeArcs);
	}
	return sc;
}

__global__ void MonteCarloArcsKernel(
    curandState * state,
    int * active_region,
    bool * clusters_fixed,
    float * clusters_positions,
    float * heatmap_exp_dist_anchor,
    int * milestone_successes,
    float * scores,
    const float T,
    const int N,
    const float score,
    const float tempJumpScale,
    const float tempJumpCoef,
    const float springConstantStretchArcs,
    const float springConstantSqueezeArcs,
    const bool use2D,
    const int numberOfParallelSimulations,
    const int step_size,
    const int numberOfClusters,
    const int activeRegionSize,
    bool * error
) {
    int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if(threadIndex >= numberOfParallelSimulations) return;

    int p, ind;
    float local_score_prev, local_score_curr;	// local score
    float tp = 0.0;							    // transition probability
    float tmp[3];
    bool ok;

    float score_curr = score;
    float score_prev = score_curr;

    curandState localState = state[threadIndex];

    for(int i = 0; i < N; ++i) {
        // p = random(size);		// select point to mutate
        p = random(activeRegionSize, &localState);       // select point to mutate
        ind = active_region[p];	                         // translate the index to 'clusters'

        if (clusters_fixed[ind]) *error = true;
        if(*error) return;

        local_score_prev = calcScoreDistancesActiveRegion(p, active_region, &(clusters_positions[threadIndex * numberOfClusters * 3]), 
                                    heatmap_exp_dist_anchor, springConstantStretchArcs, springConstantSqueezeArcs);

        randomVector(tmp, step_size, use2D, &localState);
        
        // explanation of indexing:
        // Each thread has its own dedicated state of the system, which consists of numberOfClusters elements.
        // Each element in the system is represented as 3d vector, so we need to stride by 3 array entries.
        addToVector(clusters_positions + threadIndex * numberOfClusters * 3 + ind * 3, tmp);

        local_score_curr = calcScoreDistancesActiveRegion(p, active_region, &(clusters_positions[threadIndex * numberOfClusters * 3]), 
                                    heatmap_exp_dist_anchor, springConstantStretchArcs, springConstantSqueezeArcs);

        score_curr = score_curr - local_score_prev + local_score_curr;

        ok = score_curr <= score_prev;

        if (!ok) {
            tp = tempJumpScale * expf(-tempJumpCoef * (score_curr/score_prev) / T);
            ok = withChance(tp, &localState);
        }

        if (ok) {
            milestone_successes[threadIndex]++;
        } else {
            subtractFromVector(&(clusters_positions[threadIndex * numberOfClusters * 3 + ind * 3]), tmp);
            score_curr = score_prev;
        }

        score_prev = score_curr;
    }
    state[threadIndex] = localState;
    scores[threadIndex] = score_curr;
}

__global__ void setupKernel(curandState * state, time_t seed) {
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    curand_init(seed, tid, 0, &state[tid]);
}

float LooperSolver::parallelMonteCarloArcs(float step_size) {
    const int N = 100;
    const int blocks = 32;
    const int threads = 32;
    // const int blocks = Settings::numberOfBlocks;
    // const int threads = Settings::numberOfThreads;
    const int numberOfParallelSimulations = blocks * threads;
    
	float dt = Settings::dtTemp;				// change of temperature
	float T = Settings::maxTemp;			    // temperature
    int size = active_region.size();

    int iterations = 0;

	float score_curr;		// global score (used to check if we can stop)
	float milestone_score;		        // score measured every N steps, used to see whether the solution significantly improves
    int milestone_success = 0;	        // calc how many successes there were since last milestone

	// the idea to keep the true score without recalculating it every time is to calculate total score once, on the beginning.
	// during the simulation we calculate local score before and after a move
	// new score is then equal to (old total score) - (local score before move) + (local score after move)

	// calculate total score
	score_curr = calcScoreDistancesActiveRegion();

    milestone_score = score_curr;
    size_t heatmapSize = heatmap_exp_dist_anchor.size;
    
    bool * d_hasError;
    bool h_hasError;

    cudaMalloc((void**)&d_hasError, sizeof(bool));
    cudaMemset(d_hasError, 0, sizeof(bool));

    thrust::host_vector<bool> h_clusters_fixed(clusters.size());
    thrust::host_vector<bool> h_clusters_positions(clusters.size() * 3); // just one initial state to copy across
    thrust::host_vector<float> h_scores(numberOfParallelSimulations);
    thrust::host_vector<int> h_milestone_successes(numberOfParallelSimulations);

    for(int i = 0; i < clusters.size(); ++i) {
        h_clusters_fixed[i] = clusters[i].is_fixed;

        for(int j = 0; j < 3; ++j) {
            h_clusters_positions[i * 3 + j] = clusters[i].pos[j];
        }
    }

    thrust::device_vector<int> d_active_region(active_region);
    thrust::device_vector<bool> d_clusters_fixed = h_clusters_fixed;
    thrust::device_vector<float> d_clusters_positions(numberOfParallelSimulations * clusters.size() * 3);
    thrust::device_vector<float> d_heatmap_exp_dist_anchor(size * size);
    thrust::device_vector<int> d_milestone_successes(numberOfParallelSimulations, 0); 
    thrust::device_vector<float> d_scores(numberOfParallelSimulations, 0.0);

    // TODO: make sure it's copying in the correct order
    for(int i = 0; i < heatmapSize; ++i) {
        thrust::copy(heatmap_exp_dist_anchor.v[i], heatmap_exp_dist_anchor.v[i] + heatmapSize, d_heatmap_exp_dist_anchor.begin() + heatmapSize * i);
    }
    
    for(int i = 0; i < numberOfParallelSimulations; ++i) {
        thrust::copy(h_clusters_positions.begin(), h_clusters_positions.end(), d_clusters_positions.begin() + i * clusters.size() * 3);
    }

    curandState * devStates;
	gpuErrchk(cudaMalloc((void**)&devStates, numberOfParallelSimulations * sizeof(curandState)));

    // cuRand initialization
	setupKernel<<<blocks, threads>>>(devStates, time(NULL));

    // iterators used to retrieve the best current system state
    thrust::detail::normal_iterator<thrust::device_ptr<float>> iterStart, iterEnd;

	output(4, "size = %d; initial score = %lf, step size = %f\n", size, score_curr, step_size);

    while(true) {
        // launch kernel computing N sized Markov Chain
        MonteCarloArcsKernel<<<blocks, threads>>>(
            devStates,
            thrust::raw_pointer_cast(d_active_region.data()),
            thrust::raw_pointer_cast(d_clusters_fixed.data()),
            thrust::raw_pointer_cast(d_clusters_positions.data()),
            thrust::raw_pointer_cast(d_heatmap_exp_dist_anchor.data()),
            thrust::raw_pointer_cast(d_milestone_successes.data()),
            thrust::raw_pointer_cast(d_scores.data()),
            T,
            N,
            score_curr,
            Settings::tempJumpScale,
            Settings::tempJumpCoef,
            Settings::springConstantStretchArcs,
            Settings::springConstantSqueezeArcs,
            Settings::use2D,
            numberOfParallelSimulations,
            step_size,
            clusters.size(),
            active_region.size(),
            d_hasError
        );

        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        cudaMemcpy(&h_hasError, d_hasError, sizeof(bool), cudaMemcpyDeviceToHost);
        if(h_hasError) error("cluster fixed during arcs!\n");

        int resultIndex = thrust::min_element(thrust::device, d_scores.begin(), d_scores.end()) - d_scores.begin();
        h_scores = d_scores;
        h_milestone_successes = d_milestone_successes;

        score_curr = h_scores[resultIndex];
        milestone_success = h_milestone_successes[resultIndex];

        iterStart = d_clusters_positions.begin() + resultIndex * clusters.size() * 3;
        iterEnd = d_clusters_positions.begin() + resultIndex * clusters.size() * 3 + clusters.size() * 3;

        for(int i = 0; i < numberOfParallelSimulations - 1; ++i) {
            thrust::copy(iterStart, iterEnd, d_clusters_positions.begin() + i * clusters.size() * 3);
        }

        T *= dt;  // cooling
        iterations += N;

        // TODO: make sure that Settings::MCstopConditionSteps is divisible by N
        // check if we should stop
		if (iterations % (Settings::MCstopConditionSteps / 10) == 0) {
			output(7, "milestone: score = %lf, last = %lf (%lf), T=%lf last = %d\n", score_curr, milestone_score,
					score_curr/milestone_score, T, milestone_success);
            printf("WORKS YO \n");
			// stop if the improvement since the last milestone is less than threshold,
			// or if the score is too small (may happen when heatmap is small and points are perfectly arranged)
			if ((score_curr > Settings::MCstopConditionImprovement * milestone_score &&
					milestone_success < Settings::MCstopConditionMinSuccesses) || score_curr < 1e-5 || score_curr/milestone_score>0.9999) break;

            milestone_score = score_curr;
            thrust::fill(d_milestone_successes.begin(), d_milestone_successes.end(), 0);
        }
    }
    // set best state to LooperSolver::clusters::pos
    thrust::copy(iterStart, iterEnd, h_clusters_positions.begin());

    for(int i = 0; i < clusters.size(); ++i) {
        for(int j = 0; j < 3; ++j) {
            clusters[i].pos[j] = h_clusters_positions[i * 3 + j];
            printf(" %f ", clusters[i].pos[j]);
        }
        printf("\n");
    }

    cudaFree(d_hasError);
    cudaFree(devStates);

    return score_curr;
}

__global__ void MonteCarloHeatmapKernel(
    curandState * state,
    bool * clusters_fixed,
	float * clusters_positions,
	int * heatmap_chromosome_boundaries,
	int * milestone_successes,
	float * scores,
	const float _T,
    const int N,
    const float score,
	gpu_settings settings,
	float * heatmap_dist,
    const int step_size,
    const int numberOfClusters,
    const int activeRegionSize,
    const int heatmapSize,
    const int heatmapDiagonalSize,
    const int chromosomeBoundariesSize,
    bool * error
) {
	int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if(threadIndex >= activeRegionSize) return;

    // COPY POSITIONS TO SHARED MEMORY
    extern __shared__ float d[];
    float * positions = d;
    bool * ok = (bool*) (d + activeRegionSize * 3);
    float T = _T;
    int iterations = 0;
    float milestone_score = score;
    #pragma unroll
    for(int i = 0; i < 3; ++i) {
        positions[threadIdx.x*3 + i] = clusters_positions[threadIndex * 3 + i];
    }
	float tp = 0.0;
	float displacement[3];
    bool localOk;
	float score_curr = score;
	float score_prev = score_curr;

	curandState localState = state[threadIndex];

    while(true) {
        for(int i = 0; i < N; ++i) {
            if (clusters_fixed[threadIndex]) *error = true;
            if(*error) return;
            randomVector(displacement, step_size / 15, settings.use2D, &localState);
            addToVector(positions + threadIdx.x * 3, displacement);
            __syncthreads();
            if(threadIdx.x == 0) {
                score_curr = calcScoreHeatmapActiveRegion(-1, positions, heatmap_chromosome_boundaries, 
                        heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize, chromosomeBoundariesSize);
                *ok = score_curr <= score_prev;

                if (!(*ok) && T > 0.0) {
                    tp = settings.tempJumpScaleHeatmap * expf(-settings.tempJumpCoefHeatmap * (score_curr/score_prev) / T);
                    *ok = withChance(tp, &localState);
                }

                if (*ok) {
                    milestone_successes[0]++;
                }
                else {
                    score_curr = score_prev;
                }
            }
            __syncthreads();
            if(!(*ok)) {
                subtractFromVector(&(positions[threadIndex * 3]), displacement);
            }
            score_prev = score_curr;
        }
        T *= 0.99;
        iterations += N;
        // check if we should stop
        // make sure that Settings::MCstopConditionSteps is divisible by N
        if (iterations % 200 == 0 && threadIdx.x == 0) {
            //if ((score_curr > Settings::MCstopConditionImprovementHeatmap * milestone_score &&
             //       milestone_success < Settings::MCstopConditionMinSuccessesHeatmap) || score_curr < 1e-6) {
                *ok = true;
                scores[0] = score_curr;
            //}
            //else {
            //    *ok = false;
            //    milestone_score = score_curr;
           // }
            
        }
        else if(threadIdx.x == 0) {
            *ok = false;
        }
        __syncthreads();
        if(*ok) {
            break;
        }
    }
    state[threadIndex] = localState;
    #pragma unroll
    for(int i = 0; i < 3; ++i) {
        clusters_positions[threadIndex * 3 + i] = positions[threadIndex * 3 + i];
    }
}

float LooperSolver::ParallelMonteCarloHeatmap(float step_size) {
	step_size *= 0.5;

	const int N = 100;
    const int blocks = 32;
    const int threads = 32;
    // const int blocks = Settings::numberOfBlocks;
    // const int threads = Settings::numberOfThreads;
    const int numberOfParallelSimulations = blocks * threads;

	double T = Settings::maxTempHeatmap;		// set current temperature to max

	double score_curr;
	double score_density;
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
    thrust::device_vector<float> d_clusters_positions(active_region.size() * 3);
	thrust::device_vector<int> d_milestone_successes(1, 0); 
    thrust::device_vector<float> d_scores(1, 0.0);
    thrust::device_vector<float> d_heatmap_dist(heatmapSize * heatmapSize);

    // TODO: make sure it's copying in the correct order
    for(int i = 0; i < heatmapSize; ++i) {
        thrust::copy(heatmap_dist.v[i], heatmap_dist.v[i] + heatmapSize, d_heatmap_dist.begin() + heatmapSize * i);
    }

    d_clusters_positions = h_clusters_positions;
	curandState * devStates;
	gpuErrchk(cudaMalloc((void**)&devStates, active_region.size() * sizeof(curandState)));

    // cuRand initialization
	setupKernel<<< 1, active_region.size() >>>(devStates, time(NULL));
	gpuErrchk( cudaDeviceSynchronize() );

	output(2, "initial score: %lf (density=%lf)\n", score_curr, score_density);

    MonteCarloHeatmapKernel<<<1, active_region.size(), active_region.size() * 3 * sizeof(float) + sizeof(bool)>>>(
        devStates,
        thrust::raw_pointer_cast(d_clusters_fixed.data()),
        thrust::raw_pointer_cast(d_clusters_positions.data()),
        thrust::raw_pointer_cast(d_heatmap_chromosome_boundaries.data()),
        thrust::raw_pointer_cast(d_milestone_successes.data()),
        thrust::raw_pointer_cast(d_scores.data()),
        T,
        N,
        score_curr,
        settings,
        thrust::raw_pointer_cast(d_heatmap_dist.data()),
        step_size,
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
    h_scores = d_scores;
    h_milestone_successes = d_milestone_successes;
    score_curr = h_scores[0];
    printf("FINAL SCORE IS %f", score_curr);
    milestone_success = h_milestone_successes[0];
    h_clusters_positions = d_clusters_positions;

    for(int i = 0; i < active_region.size(); ++i) {
        #pragma unroll
        for(int j =0; j < 3; ++j) {
            clusters[active_region[i]].pos[j] = h_clusters_positions[i*3 + j];
        }
    }
    printf("CPU SCORE IS %f \n", calcScoreHeatmapActiveRegion());
	// TODO: do we need more cudaFree?
	cudaFree(d_hasError);
    cudaFree(devStates);
	return score_curr;
}
