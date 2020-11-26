//TODO Move Curand to other kernel

#include <curand_kernel.h>
#include <assert.h>
#include <iostream>

struct gpu_density {
	bool* odw;	
	int size_x, size_y, size_z;
	bool is_static;
	const int dx[6] = {-1, 1, 0, 0, 0, 0};
	const int dy[6] = {0, 0, 1, -1, 0, 0};
	const int dz[6] = {0, 0, 0, 0, 1, -1};
	float* t;
	float center[3];
	float origin[3];
};

struct key_3d {
	int x, y, z;
	float value;
	key_3d* next = nullptr;
};


struct gpu_settings {
	float densityInfluence;
	float densityScale;
	float densityWeight;
	double maxTempHeatmapDensity;
	float dtTempHeatmapDensity;
	float dtTempHeatmap;
	double epsilon;
	float MCstopConditionImprovementHeatmapDensity;
	float MCstopConditionImprovementHeatmap;
	bool use2D;
	float tempJumpCoefHeatmapDensity;
	float tempJumpCoefHeatmap;
	float tempJumpScaleHeatmapDensity;
	float tempJumpScaleHeatmap;
	int MCstopConditionMinSuccessesHeatmapDensity;
	int MCstopConditionMinSuccessesHeatmap;
	int MCstopConditionStepsHeatmapDensity;
	int MCstopConditionStepsHeatmap;
	float maxTempHeatmap;
};

struct queue {
	key_3d* head = nullptr;
	key_3d* tail = nullptr;
	int size = 0;
};

__device__ void queue_push(queue q, key_3d* elem) {
	if(q.head == nullptr) {
		q.head = elem;
		q.tail = elem;
		q.size = 1;
	} else {
		q.tail->next = elem;
		q.tail = elem;
		++(q.size);
	}
}

__device__ key_3d queue_front(queue q) {
	return *(q.tail);
}

__device__ void queue_pop(queue q) {
	if(q.head == nullptr) {
		return;
	}
	else if(q.head->next == nullptr){
		q.head = nullptr;
		q.tail = nullptr;
		q.size = 0;
	} else {
		q.head = q.head->next;
		--(q.size);
	}
}
 
__device__ float vector3Length(float *vec) {
	return sqrtf(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

__device__ bool getOdw(const gpu_density &density, const int &x, const int &y, const int &z) {
	return density.odw[z * density.size_x * density.size_y + y * density.size_x + x];
}

__device__ void setOdw(const gpu_density &density, const int &x, const int &y, const int &z, const bool &val) {
	density.odw[z * density.size_x * density.size_y + y * density.size_x + x] = val;
}

__device__ float getT(const gpu_density &density, const int &x, const int &y, const int &z) {
	return density.t[z * density.size_x * density.size_y + y * density.size_x + x];
}

__device__ void setT(const gpu_density &density, const int &x, const int &y, const int &z, const float &val) {
	density.t[z * density.size_x * density.size_y + y * density.size_x + x] = val;
}

__device__ void clearOdw(gpu_density &d) {
	if (d.is_static) return;
	for (int i = 0; i < d.size_x; ++i) {
		for (int j = 0; j < d.size_y; ++j) {
			for (int k = 0; k < d.size_z; ++k) setOdw(d, i, j, k, false);
		}
	}
}

__device__ void normalize(float vec[3]) 
{
  float m = vector3Length(vec);
  if (m > 0.0F) 
    m = 1.0F / m;
  else 
    m = 0.0F;
  vec[0] *= m;
  vec[1] *= m;
}


__device__ void addPointMass(gpu_density &d, const int &x,const int &y, const int &z, const float &val, const gpu_settings &settings) {
	if (x < 0 || y < 0 || z < 0 || x >= d.size_x || y >= d.size_y || z >= d.size_z) {
		return;
	}

	queue q;
	key_3d p, np;
	key_3d first = {x, y, z, val, nullptr};
    clearOdw(d);
    // TODO FIX QUEUE PROBLEM
	queue_push(q, &first);
	setOdw(d, x, y, z, true);
	while (q.size) {
		p = queue_front(q);
		queue_pop(q);
		setT(d, p.x, p.y, p.z, max(getT(d, p.x, p.y, p.z), p.value));
		#pragma unroll
		for (int i=0; i<6; i++) {
			np.x = p.x + d.dx[i];
			np.y = p.y + d.dy[i];
			np.z = p.z + d.dz[i];
			np.value = p.value * settings.densityInfluence;
			np.next = nullptr;

			if (np.x >= 0 && np.y >= 0 && np.z >= 0 && np.x < d.size_x && np.y < d.size_y && np.z < d.size_z && np.value > 0.001f) {
				if (!getOdw(d, np.x, np.y, np.z)) {
					queue_push(q, &np);
					setOdw(d, np.x, np.y, np.z, true);
				}
			}
		}
	}
}

__device__ int random(curandState_t *s, int range) {
	return curand_uniform(s) * range;
}

__device__ float random_uniform(curandState_t *s) {
	return curand_uniform(s);
}

__device__ float random(curandState_t *s, float range, bool negative) {
	if (negative) return (2.0f * random_uniform(s) - 1.0f) * range;
	return range * random_uniform(s);
}

__device__ void random_vector(float vec[3], curandState_t *s, float max_size, bool in2D = false) {
	if (in2D) {
		#pragma unroll
		for(int i = 0; i < 2; ++i)
			vec[i] = random(s, max_size, true);
		vec[2] = 0.0;
	}
	else {
		#pragma unroll
		for(int i = 0; i < 3; ++i)
				vec[i] = random(s, max_size, true);
	}
}

__device__ bool withChance(curandState_t *s, float chance) {
	return random_uniform(s) < chance;
}

// should be okay
__device__ void normalize(gpu_density &d, const gpu_settings &settings) {
	float total = 0.0f;
	for (int i = 0; i < d.size_x; ++i) {
		for (int j = 0; j < d.size_y; ++j) {
			for (int k = 0; k < d.size_z; ++k) total += getT(d, i, j, k);
		}
	}

	float avg = total / ((float)d.size_x * d.size_y * d.size_z);
	if (avg < settings.epsilon) return;
	for (int i = 0; i < d.size_x; ++i) {
		for (int j = 0; j < d.size_y; ++j) {
			for (int k = 0; k < d.size_z; ++k) setT(d, i, j, k, getT(d, i, j, k) / avg);
		}
	}
}

__device__ void addToVector3(float *destination, const float *value1, const float *value2) {
	#pragma unroll
    for(int i = 0; i < 3; ++i) {
        destination[i] = value1[i] - value2[i];
    }
}

__device__ void subtractFromVector3(float *destination, const float *value1, const float *value2) {
	#pragma unroll
    for(int i = 0; i < 3; ++i) {
        destination[i] = value1[i] - value2[i];
    }
}

__device__ void scaleVector3(float *destination, const float &value) {
	#pragma unroll
	for(int i = 0; i < 3; ++i) {
		destination[i] *= value;
	}
}



__device__ void clear_density(gpu_density &density_curr) {
	for (int i = 0; i < density_curr.size_x; ++i) {
		for (int j = 0; j < density_curr.size_y; ++j) {
			for (int k = 0; k < density_curr.size_z; ++k) {
				setT(density_curr, i, j, k, 0.0f);
				if (!density_curr.is_static) setOdw(density_curr, i, j, k, false);
			}
		}
	}
}

__device__ void repairDensityBoundary(gpu_density &density_curr, int* active_region, gpu_density &density, float* clusters_positions, const gpu_settings &settings, curandState_t* cur_state) {
	clear_density(density_curr);
	size_t n = sizeof(active_region) / sizeof(active_region[0]);
	float pos[3];
	float shift[3];
	subtractFromVector3(shift, density.center, density.origin);
	int px, py, pz;
	for (size_t i = 0; i < n; ++i) {
		while (1) {
			subtractFromVector3(pos, clusters_positions + (active_region[i] * 3), density.center);
			scaleVector3(pos, settings.densityScale);
			pos[0] /= 3.0f;
			addToVector3(pos, pos, shift);
			px = (int)pos[0];
			py = (int)pos[1];
			pz = (int)pos[2];
			if (px < 0 || py < 0 || pz < 0 || px >= density.size_x || py >= density.size_y || pz >= density.size_z || getT(density, px, py, pz) < settings.epsilon) {
				float cen[3];
				float ran[3];
				random_vector(ran, cur_state, 20.0f);
				addToVector3(cen, density.center, ran);
				float shift[3];
				subtractFromVector3(shift, cen, clusters_positions + (active_region[i] * 3));
				normalize(shift);
				scaleVector3(shift, 2.0f);
				addToVector3(clusters_positions + (active_region[i] * 3), clusters_positions + (active_region[i] * 3), shift);
			}
			else break;
		}
	}

	normalize(density_curr, settings);
}


//this function uses following external objects and functionalities:
// -> density (uses non changing fields - we can easily just copy that to kernel memory and use the same on all threads)
//		-> center
//		-> origin
//		-> size_x
//		-> size_y
//		-> size_z
//		-> 3 dimensional vector of floats t, each thread is going to use a all elements but in read - only mode
// -> density functions such as: 
//		-> addPointMass	
//		-> normalize
// -> size of the active region (read only)
// -> the whole active region
// -> pos field of active region index in the clusters array (read only)
// -> epsilon
// VECTOR3
//	->.x .y .z
//	-> length function
// overloaded addition
// settings object
__device__ double calcScoreDensity(const gpu_density &density, gpu_density &density_curr, int* active_region, const float* clusters_positions, const gpu_settings &settings) {
	double error = 0.0, e;
	clear_density(density_curr);
	size_t n = sizeof(active_region)/sizeof(active_region[0]);;
	float pos[3];
	float shift[3];
	subtractFromVector3(shift, density.center, density.origin);
	int px, py, pz;
	for (size_t i = 0; i < n; ++i) {
        // find 3d density-coordinates for the current bead
        const float * local_ptr = clusters_positions + (active_region[i] * 3);
		subtractFromVector3(pos, local_ptr, density.center);
		scaleVector3(pos, settings.densityScale);
		pos[0] /= 3.0f;
		addToVector3(pos, pos, shift);
		px = (int)pos[0];
		py = (int)pos[1];
		pz = (int)pos[2];

		if (px < 0 || py < 0 || pz < 0 || px >= density.size_x || py >= density.size_y || pz >= density.size_z) {
			int d = 0;
			if (px < 0) d += -px;
			else if (px >= density.size_x) d += px - density.size_x + 1;
			if (py < 0) d += -py;
			else if (py >= density.size_y) d += py - density.size_y + 1;
			if (pz < 0) d += -pz;
			else if (pz >= density.size_z) d += pz - density.size_z + 1;
			error += d * 1000000;
		}
		else if (getT(density, px, py, pz) < settings.epsilon) {
            float temp[3];
            subtractFromVector3(temp, pos, density.center);
			error += 1000000.0f +  1000000.0f / vector3Length(temp);
		}
		else {
			addPointMass(density_curr, px, py, pz, 1.0f, settings);
		}
	}

	normalize(density_curr, settings);

	for (int i = 0; i < density.size_x; ++i) {
		for (int j = 0; j < density.size_y; ++j) {
			for (int k = 0; k < density.size_z; ++k) {
				if (getT(density, i, j, k) > settings.epsilon) {
					e = (getT(density, i, j, k) - getT(density_curr, i, j, k)) / getT(density, i, j, k);
					e = 0.0;
				}
				else {
					e = 5.0f * getT(density_curr, i, j, k);
				}
				error += e * e;
			}
		}
	}
	return error * settings.densityWeight;
}

//this function uses following external funcitonality / data
// -> size of heatmap_chromosome_boundaries
// -> the entire heatmap_chromosome_boundaries
__device__ void getChromosomeHeatmapBoundary(int p, int &start, int &end, int* gpu_heatmap_chromosome_boundaries) {
	size_t hcb_size = sizeof(gpu_heatmap_chromosome_boundaries) / sizeof(gpu_heatmap_chromosome_boundaries[0]);
	if (hcb_size == 0) return;
	if (p < 0 || p > gpu_heatmap_chromosome_boundaries[hcb_size - 1]) return;

	for (size_t i = 0; i+1 < hcb_size; ++i) {
		if (gpu_heatmap_chromosome_boundaries[i] <= p && p < gpu_heatmap_chromosome_boundaries[i+1]) {
			start = gpu_heatmap_chromosome_boundaries[i];
			end = gpu_heatmap_chromosome_boundaries[i+1]-1;
			return;
		}
	}
}

struct gpu_heatmap_dist {
	int diagonal_size;
    float** v;
};


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
__device__ double calcScoreHeatmapActiveRegion(int moved, float* clusters_positions, int* active_region, int* gpu_heatmap_chromosome_boundaries, const gpu_heatmap_dist &heatmap_dist) {
	double err = 0.0, cerr;
	float temp[3];
	float d;
	size_t n = sizeof(active_region)/sizeof(active_region[0]);;
	size_t hd_size = sizeof(heatmap_dist.v) / sizeof(heatmap_dist.v[0]);
	size_t hcb_size = sizeof(gpu_heatmap_chromosome_boundaries) / sizeof(gpu_heatmap_chromosome_boundaries[0]);
	if (hd_size != n) { printf("heatmap sizes mismatch, dist size=%zu, active region=%zu", hd_size, n); return 0.0; }
	if (moved == -1) {
		for (size_t i = 0; i < n; ++i) err += calcScoreHeatmapActiveRegion(i, clusters_positions, active_region, gpu_heatmap_chromosome_boundaries, heatmap_dist);
	}
	else {
		int st = 0;
		int end = n- 1;
		if (hcb_size > 0) getChromosomeHeatmapBoundary(moved, st, end, gpu_heatmap_chromosome_boundaries);

		for (int i = st; i <= end; ++i) {
			if (abs(i-moved) >= heatmap_dist.diagonal_size) {
                if (heatmap_dist.v[i][moved] < 1e-6) continue;	// ignore 0 values
                const float* temp_one = clusters_positions + active_region[i] * 3;
                const float* temp_two = clusters_positions + active_region[moved] * 3;
				subtractFromVector3(temp, temp_one, temp_two);
				d = vector3Length(temp);
				cerr = (d - heatmap_dist.v[i][moved]) / heatmap_dist.v[i][moved];
				err += cerr * cerr;
			}
		}
	}
	return err;
}


__device__ void initialize_density(gpu_density &density_new, const gpu_density &density_old, float* t_space, bool * odw_space, int idx) {
	density_new.odw = odw_space + density_old.size_x * density_old.size_y * density_old.size_z * idx;	
	density_new.size_x = density_old.size_x;
	density_new.size_y = density_old.size_y;
	density_new.size_z = density_old.size_z;
	density_new.is_static = false;
    density_new.t = t_space + density_old.size_x * density_old.size_y * density_old.size_z * idx;
    #pragma unroll
    for(int i = 0; i < 3; ++i) {
        density_new.center[i] = density_old.center[i];
	    density_new.origin[i] = density_old.origin[i];
    }
	
}



// this function operates on data such as:
// -> active_region
//		-> size
//	VECTOR3 structure and functions as described above
// the whole random issue
// entire active region
// clusters
//		->is_fixed
//		->pos
// clusters[i].pos is edited and adjusted by thread to get better state, we will have to exchange it
// settings object
// withChance method from CPP (better use curand)

__global__ void MonteCarloHeatmapAndDensity(float step_size,
											int* active_region, 
											bool* clusters_is_fixed, 
											float* clusters_positions, 
											int state_size,
											gpu_settings settings, 
											gpu_heatmap_dist heatmap_dist,
											gpu_density density,
											bool* density_bool,
                                            float* density_float,
                                            int* gpu_heatmap_chromosome_boundaries) {
    int idx_g = threadIdx.x + blockDim.x * blockIdx.x;
	float* local_clusters_positions = clusters_positions + idx_g * 3 * state_size;
	double T = settings.maxTempHeatmapDensity;		// set current temperature to max
	double tp = 0.0;	// transition probability
	float displacement[3];
	int i, p, ind;
	bool ok;
	double score_prev, score_curr;		// global score (used to check if we can stop)
	double score_density, score_heatmap;
	double score_density_prev, score_heatmap_prev;
	double milestone_score;		// score measured every N steps, used to see whether the solution significantly improves
	int success = 0;			// overall number of successes
	int milestone_success = 0;	// calc how many successes there were since last milestone
	int size = sizeof(active_region)/sizeof(active_region[0]);
	if (size <= 1) return;	// there is nothing to do
	gpu_density density_curr;
	initialize_density(density_curr, density, density_float, density_bool, idx_g);
	//initialize curand
	curandState_t* cur_state;
	curand_init(2019UL, idx_g, 0, cur_state);

	// calc initial score
	score_heatmap = calcScoreHeatmapActiveRegion(-1, local_clusters_positions, active_region, gpu_heatmap_chromosome_boundaries, heatmap_dist);	// score from heatmap
	score_density = calcScoreDensity(density, density_curr, active_region, local_clusters_positions, settings);
	score_curr = score_heatmap + score_density;
	//score_curr = score_density;
	score_prev = score_curr;
	milestone_score = score_curr;
	int milestone_cnt = 0;
	i = 1;
	while (true) {
		// select point to mutate and find chromosome on which it is located
		p = random(cur_state, size);	// index in 'active_region'
		ind = active_region[p];	// index in 'clusters'
		if (clusters_is_fixed[ind]) continue;	// check if bead is fixed (for example telomere)
		random_vector(displacement, cur_state, step_size, settings.use2D);	// generate random displacement vector
		const float* temp_one = local_clusters_positions + ind * 3;
		addToVector3(local_clusters_positions + ind * 3, temp_one , displacement);
		score_heatmap = calcScoreHeatmapActiveRegion(p, local_clusters_positions, active_region, gpu_heatmap_chromosome_boundaries, heatmap_dist);
		score_density = calcScoreDensity(density, density_curr, active_region, local_clusters_positions, settings);
		score_curr = score_heatmap + score_density;
		ok = score_curr <= score_prev;
		if (!ok && T > 0.0) {
			tp = settings.tempJumpScaleHeatmapDensity * expf(-settings.tempJumpCoefHeatmapDensity * (score_curr/score_prev) / T);
			ok = withChance(cur_state, tp);
		}
		if (ok) {
			success++;
			milestone_success++;
		}
		else {
			subtractFromVector3(local_clusters_positions + ind * 3, temp_one, displacement);
			score_curr = score_prev;		// score doesn't change
			score_density = score_density_prev;
			score_heatmap = score_heatmap_prev;
		}
		T *= settings.dtTempHeatmapDensity;    // cooling
		// check if we should stop
		if (i % settings.MCstopConditionStepsHeatmapDensity == 0) {
			// stop if the improvement since the last milestone is less than 0.5%, or if the score is too small (may happen when heatmap is small and points are perfectly arranged)
			if ((score_curr > settings.MCstopConditionImprovementHeatmapDensity * milestone_score &&
					milestone_success < settings.MCstopConditionMinSuccessesHeatmapDensity) || score_curr < 1e-6) {
				break;
			}
			milestone_score = score_curr;
			milestone_success = 0;
			milestone_cnt++;
			repairDensityBoundary(density_curr, active_region, density, local_clusters_positions, settings, cur_state); 
		}
		score_prev = score_curr;
		score_heatmap_prev = score_heatmap;
		score_density_prev = score_density;
		++i;
	}
	// return score_curr;
}

__global__ void MonteCarloHeatmap(  float step_size,
									int* active_region, 
									bool* clusters_is_fixed, 
									float* clusters_positions, 
									int state_size,
									gpu_settings settings, 
									gpu_heatmap_dist heatmap_dist,
									gpu_density density,
									bool* density_bool,
									float* density_float,
									int *gpu_heatmap_chromosome_boundaries) {

	step_size *= 0.5;
	double T = settings.maxTempHeatmap;		// set current temperature to max
    int idx_g = threadIdx.x + blockDim.x * blockIdx.x;
	double tp = 0.0;	// transition probability
	float displacement[3];
	int i, p, ind;
	bool ok;
	double score_prev, score_curr;		// global score (used to check if we can stop)
	//TODO Density?
	//double score_density;
	double milestone_score;		// score measured every N steps, used to see whether the solution significantly improves
	int success = 0;			// overall number of successes
	int milestone_success = 0;	// calc how many successes there were since last milestone
	int size = sizeof(active_region) / sizeof(active_region[0]);
	if (size <= 1) return;	// there is nothing to do return 0
	float* local_clusters_positions = clusters_positions + idx_g * 3 * state_size;
	gpu_density density_curr;
	initialize_density(density_curr, density, density_float, density_bool, idx_g);
	// calc initial score
	score_curr = calcScoreHeatmapActiveRegion(-1, local_clusters_positions, active_region, gpu_heatmap_chromosome_boundaries, heatmap_dist);	// score from heatmap
	//TODO decide what's going to happen to this one
	//score_density = calcScoreDensity(density, density_curr, active_region, local_clusters_positions, settings);
	score_prev = score_curr;
	milestone_score = score_curr;
	int milestone_cnt = 0;
	curandState_t* cur_state;
	curand_init(2019UL, idx_g, 0, cur_state);
	i = 1;
	while (true) {
		// select point to mutate and find chromosome on which it is located
		p = random(cur_state, size);	// index in 'active_region'
		ind = active_region[p];	// index in 'clusters'
		if (clusters_is_fixed[ind]) continue;	// check if bead is fixed (for example telomere)
		random_vector(displacement, cur_state, step_size, settings.use2D);	// generate random displacement vector
		addToVector3(local_clusters_positions + ind * 3, local_clusters_positions + ind * 3, displacement);
		score_curr = calcScoreHeatmapActiveRegion(p, local_clusters_positions, active_region, gpu_heatmap_chromosome_boundaries, heatmap_dist);
		ok = score_curr <= score_prev;
		if (!ok && T > 0.0) {
			tp = settings.tempJumpScaleHeatmap * expf(-settings.tempJumpCoefHeatmap * (score_curr/score_prev) / T);
			ok = withChance(cur_state, tp);
		}
		if (ok) {
			success++;
			milestone_success++;
		}
		else {
			subtractFromVector3(local_clusters_positions + ind * 3, local_clusters_positions + ind * 3, displacement);
			score_curr = score_prev;		// score doesn't change
		}
		T *= settings.dtTempHeatmap;    // cooling
		// check if we should stop
		if (i % settings.MCstopConditionStepsHeatmap == 0) {
			//TODO another density to decide whether stays or goes away
			//score_density = calcScoreDensity(density, density_curr, active_region, local_clusters_positions, settings);
			// stop if the improvement since the last milestone is less than 0.5%, or if the score is too small (may happen when heatmap is small and points are perfectly arranged)
			if ((score_curr > settings.MCstopConditionImprovementHeatmap * milestone_score &&
					milestone_success < settings.MCstopConditionMinSuccessesHeatmap) || score_curr < 1e-6) {
				break;
			}
			milestone_score = score_curr;
			milestone_success = 0;
			milestone_cnt++;
		}
		score_prev = score_curr;
		i++;
	}
	//TODO save score_curr and score_density
}

int main() {
    printf("ALL GOOOD BROTHER\n");
    return 0;
}