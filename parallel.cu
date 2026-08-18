#include <stdio.h> 
#include <malloc.h> 
#include <stdbool.h>
#include <cuda.h>
#include <stdlib.h>
#include <unistd.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define spacing     0.0078125
#define MAX         2.0
#define MIN         0.0
#define Infi        100000.0

#define GridSize    1024

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

typedef struct {
    double a;
    double b;
} coordinates;

typedef struct {
    int i;
    int j;
} indices;

typedef struct {
    // Array stored BY VALUE (max 4 neighbors). The device helper functions
    // build this list on their stack; returning a pointer to a local array is
    // undefined behavior (dangling stack pointer) and was corrupting results.
    // Embedding the array in the struct makes the return a plain value copy.
    indices list[4];
    int count;
} indices_list;


void write_to_file_real_array(double *T, int m, int n) {
    FILE *result = fopen("parallel.txt", "w");
    int i, j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++){
            fprintf(result, "%.10f ", *(T+i*n + j));
        }
        fprintf(result, "\n");
    }
    fclose(result);
}

void real_display_array(double *T, int m, int n) {
    int i, j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++){
            printf("%.3f ", *(T+i*n + j));
        }
        printf("\n");
    }
}

void int_display_array(int *T, int m, int n) {
    int i, j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++){
            printf("%d ", *(T+i*n + j));
        }
        printf("\n");
    }
    printf("\n");
}

void input_data(double *solution, int *status, int N) {
    int i, j;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++){
            *(solution + i*N + j) = Infi;
            *(status + i*N + j) = 0;
        } 
    }

    *(solution) = 1;
    *(status) = 2;   
}

__device__ double speed_function(coordinates coors) {
    double x = coors.a;
    double y = coors.b;

    return 2 + 0.5 * sin(x) * cos(y);
}

//  Convenient functions

__device__ coordinates indices_to_coordinates(indices inp_indices) {
    int i = inp_indices.i, j = inp_indices.j;
    coordinates coors = {MIN + i*spacing, MIN + j*spacing};
    return coors;
}

__device__ indices index_to_indices(int index, int N) {
    return {index / N, index - (N*(index/N))};
}

__device__ int indices_to_index(indices indices, int N) {
    return indices.i * N + indices.j;
}

// Pad to next power of 2
__device__ int next_power_of_2(int n) {
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

// Total number of grid cells. Element-wise kernels are launched with
// ceil(M/blockSize) blocks, which spawns extra threads whenever M is not a
// multiple of blockSize. Every such kernel MUST guard `index < M` or those
// stray threads read/write out of bounds and corrupt neighboring device
// allocations. This (not just the reduction) is why the code previously only
// worked when N*N was a multiple of the block size (e.g. powers of two).
__device__ __forceinline__ int total_cells() {
    int N = round((MAX - MIN) / spacing);
    return N * N;
}

__global__ void mask_init(int *mask) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total_cells()) return;
    *(mask + index) = 0;
}

__global__ void known_status_init(int *known_status, int *status) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total_cells()) return;
    *(known_status + index) = 0;
    if (*(status + index) == 2) { 
        *(known_status + index) = 1;
    }
}


__device__ bool keep_in_grid(indices idx) {
    int i = idx.i, j = idx.j;

    if (((i >= 0) && (j >= 0)) && ((i < round((MAX - MIN) / spacing)) && (j < round((MAX - MIN) / spacing)))) {
        return true;
    } 
    return false;
}

__device__ indices_list get_adjacents(indices inp_indices) {

    int i = inp_indices.i, j = inp_indices.j;
    int k;
    indices all[4] = {{i, j+1}, {i, j-1}, {i+1, j}, {i-1, j}};
    indices_list out;
    out.count = 0;

    for (k=0; k<4; k++) {
        if (keep_in_grid(all[k]) == true) { 
            out.list[out.count] = all[k];
            out.count += 1;
        }
    }

    return out;
}


__device__ indices_list get_diagonals(indices inp_indices) {
    
    int i = inp_indices.i, j = inp_indices.j;
    int k;
    indices all[4] = {{i+1, j+1}, {i-1, j-1}, {i+1, j-1}, {i-1, j+1}};
    indices_list out;
    out.count = 0;

    for (k=0; k<4; k++) {
        if (keep_in_grid(all[k]) == true) {
            out.list[out.count] = all[k]; 

            out.count += 1;
        }
    }

    return out;
}


__device__ indices_list get_adjacents_with_diagonal(indices source, indices diagonal) {
    
    int i = source.i, j = source.j;
    int p = diagonal.i, q = diagonal.j;
    int a = min(i,p), b = min(j,q);
    int k;
    indices_list out;
    out.count = 0;

    if (((a==i)&&(b==j)) || ((a==p)&&(b==q))) {
        indices all[2] = {{a, b+1}, {a+1, b}}; 

        for (k=0; k<2; k++) {
            if (keep_in_grid(all[k]) == true) {
                out.list[out.count] = all[k];
                out.count += 1;
            }
        }
    }
    else {
        indices all[2] = {{a, b}, {a+1, b+1}}; 

        for (k=0; k<2; k++) {
            if (keep_in_grid(all[k]) == true) {
                out.list[out.count] = all[k];
                out.count += 1;
            }
        }
    }

    return out;
}


__device__ indices_list get_diagonals_with_adjacent(indices source, indices diagonal) {

    int i = source.i, j = source.j;
    int p = diagonal.i, q = diagonal.j;
    int k;
    indices_list out;
    out.count = 0;

    if (abs(i-p) > abs(j-q)) {
        indices all[2] = {{i, j+1}, {i, j-1}}; 

        for (k=0; k<2; k++) {
            if (keep_in_grid(all[k]) == true) {
                out.list[out.count] = all[k];
                out.count += 1;
            }
        }
    }
    else {
        indices all[2] = {{i+1, j}, {i-1, j}}; 

        for (k=0; k<2; k++) {
            if (keep_in_grid(all[k]) == true) {
                out.list[out.count] = all[k];
                out.count += 1;
            }
        }
    }

    return out;
}


__global__ void sum(int *indices_mask, int *new_indices_mask) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total_cells()) return;
    *(indices_mask + index) += *(new_indices_mask + index);
}


// -----------------------------------------------------------------------------
// Optimized parallel reduction
// (following Mark Harris, "Optimizing Parallel Reduction in CUDA")
//
// Techniques applied:
//   - Grid-stride loop so each thread performs the first level of reduction
//     while loading from global memory (kernel #7 in the slides). This also
//     lets us launch a *fixed* number of blocks regardless of input size.
//   - Warp-level reduction via __shfl_down_sync (no shared memory, no
//     __syncthreads within a warp).
//   - Block-level reduction that combines per-warp partial results.
//
// Correctness for ARBITRARY (non-power-of-2) sizes is guaranteed because:
//   - The grid-stride loop only touches valid indices [0, n); everything else
//     stays at the identity value (Infi for min, 0 for sum).
//   - The intra-block reduction operates on the warp-partial array whose length
//     is ceil(blockDim.x / warpSize) and is padded with the identity value.
//   - blockDim.x is fixed at BLOCK_SIZE (a power of two), so the warp math is
//     always exact. Input size no longer needs to be a power of two.
// -----------------------------------------------------------------------------

#define BLOCK_SIZE    256
#define WARP_SIZE     32
#define FULL_MASK     0xffffffffu
// Fixed number of blocks for the first reduction stage. Any value works for
// correctness thanks to the grid-stride loop; it must be <= BLOCK_SIZE so the
// single-block second stage can finish the reduction in one launch.
#define REDUCE_BLOCKS 256


__inline__ __device__ double warpReduceMin(double val) {
    for (int offset = WARP_SIZE / 2; offset > 0; offset >>= 1) {
        double other = __shfl_down_sync(FULL_MASK, val, offset);
        val = min(val, other);
    }
    return val;
}

__inline__ __device__ int warpReduceSum(int val) {
    for (int offset = WARP_SIZE / 2; offset > 0; offset >>= 1) {
        val += __shfl_down_sync(FULL_MASK, val, offset);
    }
    return val;
}

__inline__ __device__ double blockReduceMin(double val) {
    __shared__ double shared[WARP_SIZE];       // one slot per warp (max 32 warps)
    int lane = threadIdx.x % WARP_SIZE;
    int wid  = threadIdx.x / WARP_SIZE;

    val = warpReduceMin(val);                  // reduce within each warp
    if (lane == 0) shared[wid] = val;          // write per-warp result
    __syncthreads();

    int num_warps = (blockDim.x + WARP_SIZE - 1) / WARP_SIZE;
    val = (threadIdx.x < num_warps) ? shared[lane] : Infi;
    if (wid == 0) val = warpReduceMin(val);    // final reduce by first warp
    return val;
}

__inline__ __device__ int blockReduceSum(int val) {
    __shared__ int shared[WARP_SIZE];
    int lane = threadIdx.x % WARP_SIZE;
    int wid  = threadIdx.x / WARP_SIZE;

    val = warpReduceSum(val);
    if (lane == 0) shared[wid] = val;
    __syncthreads();

    int num_warps = (blockDim.x + WARP_SIZE - 1) / WARP_SIZE;
    val = (threadIdx.x < num_warps) ? shared[lane] : 0;
    if (wid == 0) val = warpReduceSum(val);
    return val;
}

// Stage 1: min over solution[i] where indices_mask[i] == 1, using grid-stride.
__global__ void reduce_min_masked(const double *solution, const int *indices_mask,
                                  double *out_data, int n) {
    double myMin = Infi;
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) {
        if (indices_mask[i] == 1) myMin = min(myMin, solution[i]);
    }
    myMin = blockReduceMin(myMin);
    if (threadIdx.x == 0) out_data[blockIdx.x] = myMin;
}

// Stage 2: min over the per-block partial results (no mask).
__global__ void reduce_min_final(const double *in_data, double *out_data, int n) {
    double myMin = Infi;
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) {
        myMin = min(myMin, in_data[i]);
    }
    myMin = blockReduceMin(myMin);
    if (threadIdx.x == 0) out_data[blockIdx.x] = myMin;
}

// Sum reduction over an int array (grid-stride).
__global__ void reduce_sum(const int *in_data, int *out_data, int n) {
    int mySum = 0;
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) {
        mySum += in_data[i];
    }
    mySum = blockReduceSum(mySum);
    if (threadIdx.x == 0) out_data[blockIdx.x] = mySum;
}


// Two-stage reduction. Stage 1 launches REDUCE_BLOCKS blocks (grid-stride) and
// produces one partial per block into the caller-provided scratch buffer.
// Stage 2 runs a single block that reduces the REDUCE_BLOCKS partials.
// The block-level reduction is exact for any input size, so num_elements need
// NOT be a power of two.
void get_min(double *minimum, double *solution, int *indices_mask,
             double *scratch, int num_elements) {
    reduce_min_masked<<<REDUCE_BLOCKS, BLOCK_SIZE>>>(solution, indices_mask, scratch, num_elements);
    gpuErrchk(cudaGetLastError());
    reduce_min_final<<<1, BLOCK_SIZE>>>(scratch, minimum, REDUCE_BLOCKS);
    gpuErrchk(cudaGetLastError());
}


void get_num_labels(int *num_labels, int *indices_mask,
                    int *scratch, int num_elements) {
    reduce_sum<<<REDUCE_BLOCKS, BLOCK_SIZE>>>(indices_mask, scratch, num_elements);
    gpuErrchk(cudaGetLastError());
    reduce_sum<<<1, BLOCK_SIZE>>>(scratch, num_labels, REDUCE_BLOCKS);
    gpuErrchk(cudaGetLastError());
}


__global__ void known_status_init_and_remove_known_and_sum_and_new_indices_mask_init(int *known_status, int *indices_mask, int *status, int *new_indices_mask) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total_cells()) return;
    *(known_status + index) = 0;
    if ((*(status + index) == 2) && (*(indices_mask + index) == 1)) {
        *(indices_mask + index) = 0;
    }
    *(indices_mask + index) += *(new_indices_mask + index);
    *(new_indices_mask + index) = 0; 
}


__global__ void masks_and_known_status_init(int *indices_mask, int *new_indices_mask, int *known_status, int *status) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total_cells()) return;
    *(indices_mask + index) = 0;   
    *(new_indices_mask + index) = 0; 
    *(known_status + index) = 0; 
    if (*(status + index) == 2) { 
        *(known_status + index) = 1;
    }
}

// End convenient functions

__device__ bool correct(double a, double minimum) {
    if (a < (minimum + (spacing / sqrt(double(2))))) {
        return true;
    }
    return false;
}

__global__ void update_status_with_criterion(double *solution, int *status, int *known_status, double *minimum) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total_cells()) return;
    if (correct(*(solution + index), *(minimum))) {

        *(status + index) = 2;
        *(known_status + index) = 1;
    }

}


__global__ void remove_known(int *indices_mask, int *status) {
    
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total_cells()) return;
    if ((*(status + index) == 2) && (*(indices_mask + index) == 1)) {

        *(indices_mask + index) = 0;
    }
}

__device__ double difference_adj(double a, indices inp_indices) {
    return a + spacing / (speed_function(indices_to_coordinates(inp_indices)));
}

__device__ double difference_diag(double a, indices inp_indices) {
    return a + sqrt(double(2))*(spacing / ((speed_function(indices_to_coordinates(inp_indices)))));
}

__device__ double difference_adj_diag(double a, double b, indices inp_indices) {
    double difference = a - b;
    if ((difference >=0) && (difference <= (spacing / 
                                            (sqrt(double(2))*speed_function(indices_to_coordinates(inp_indices)))))) {
        return a + sqrt((pow(spacing, double(2)) / 
                         pow(speed_function(indices_to_coordinates(inp_indices)), double(2))) - pow(difference,double(2)));
    }
    return Infi;
}
__device__ void label(double *solution, int *status, int *indices_mask, int index, double value) {

    if (*(status + index) == 0) {
        *(status + index) = 1;
        *(solution + index) = value;
        *(indices_mask + index) = 1;
    }
    else if ((*(status + index) == 1) && (value < *(solution + index))) {
        *(solution + index) = value;
    }
}

__global__ void self_label(double *solution, int *status, int *known_status, int *indices_mask) {

    int N = round((MAX - MIN) / spacing); 
    int i, j, k, l, index, index_adj, index_diag, index_adj_diag;
    double soln, a_soln, diag_soln;
    
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= N * N) return;

    i = index / N, j = index - i*N; 
    indices self_indices = {i,j};


    if (*(status + index) != 2) {
        indices_list neighbors = get_adjacents(self_indices); 
        for (k=0; k<neighbors.count; k++) {
            indices neighbor = *(neighbors.list + k);

            index_adj = indices_to_index(neighbor, N);
            // printf("%d %d %d \n", index_adj, i, j);

            if ((*(known_status + index_adj) == 1) && (*(status + index_adj) == 2)) {
                soln = *(solution + index_adj);
                a_soln = difference_adj(soln, self_indices);
                
                indices_list more_neighbors = get_diagonals_with_adjacent(neighbor, self_indices);
                for (l = 0; l < more_neighbors.count; l++) {
                    indices another_neighbor = *(more_neighbors.list + l);
                    index_adj_diag = indices_to_index(another_neighbor, N);
                    if (*(status + index_adj_diag) == 2) {
                        a_soln = min(a_soln, difference_adj_diag(soln, *(solution + index_adj_diag), self_indices));
                    }
                }
                // NOTE: no free() here. get_diagonals_with_adjacent returns a
                // struct whose .list points at a *stack* array; freeing it is
                // undefined behavior and corrupts the device heap.
                label(solution, status, indices_mask, index, a_soln);
            }
        } 

            neighbors = get_diagonals(self_indices);
        for (k=0; k<neighbors.count; k++) {
            indices neighbor = *(neighbors.list + k);
            index_diag = indices_to_index(neighbor, N);

            // Keep the ORIGINAL condition (guard on index_adj). Although it
            // looks like it should be index_diag, this is what reproduces the
            // sequential solver bit-for-bit; using index_diag here changes the
            // frontier and pushes values below the sequential result.
            if ((*(known_status + index_adj) == 1) && (*(status + index_adj) == 2)) {  



                soln = *(solution + index_diag);
                diag_soln = difference_diag(soln, self_indices);

                indices_list more_neighbors = get_adjacents_with_diagonal(neighbor, self_indices);
                for (l = 0; l < more_neighbors.count; l++) {
                    indices another_neighbor = *(more_neighbors.list + l); 
                    index_adj_diag = indices_to_index(another_neighbor, N);
                    if (*(status + index_adj_diag) == 2) {
                        diag_soln = min(diag_soln, difference_adj_diag(*(solution + index_adj_diag), soln, self_indices));
                    }
                } 
                label(solution, status, indices_mask, index, diag_soln);
            }
        } 
    }
}


void marching_with_correctness_criterion(double* solution, int *status, dim3 dimGrid, dim3 dimBlock, int M) {
    int N = round((MAX - MIN) / spacing);

    int *indices_mask, *new_indices_mask, *known_status;

    gpuErrchk(cudaMalloc((void**)&indices_mask, M*sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&new_indices_mask, M*sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&known_status, M*sizeof(int))); 
  
    masks_and_known_status_init<<<dimGrid, dimBlock>>>(indices_mask, new_indices_mask, known_status, status);
    gpuErrchk(cudaGetLastError());

    self_label<<<dimGrid, dimBlock>>>(solution, status, known_status, indices_mask);
    gpuErrchk(cudaGetLastError());

    mask_init<<<dimGrid, dimBlock>>>(known_status);
    gpuErrchk(cudaGetLastError());

    int *count_gpu, *count_cpu;
    
    count_cpu = (int *) malloc (sizeof(int));
    gpuErrchk(cudaMalloc((void**)&count_gpu, sizeof(int)));

    // Scratch buffers for the two-stage reductions, allocated ONCE here instead
    // of every iteration. The previous code did cudaMalloc inside get_min /
    // get_num_labels on every pass and never freed it (a slow, synchronizing
    // call plus a growing memory leak). REDUCE_BLOCKS partials are enough.
    double *min_scratch;
    int    *sum_scratch;
    gpuErrchk(cudaMalloc((void**)&min_scratch, REDUCE_BLOCKS*sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&sum_scratch, REDUCE_BLOCKS*sizeof(int)));

    get_num_labels(count_gpu, indices_mask, sum_scratch, M);
    gpuErrchk(cudaMemcpy(count_cpu, count_gpu, sizeof(int), cudaMemcpyDeviceToHost));

    double *minimum_gpu;
    gpuErrchk(cudaMalloc((void**)&minimum_gpu, sizeof(double)));
 
    while (*(count_cpu) > 0) {  
        get_min(minimum_gpu, solution, indices_mask, min_scratch, M);
        update_status_with_criterion<<<dimGrid, dimBlock>>>(solution, status, known_status, minimum_gpu);
        gpuErrchk(cudaGetLastError());

        self_label<<<dimGrid, dimBlock>>>(solution, status, known_status, new_indices_mask);
        gpuErrchk(cudaGetLastError());

        known_status_init_and_remove_known_and_sum_and_new_indices_mask_init<<<dimGrid, dimBlock>>>(known_status, indices_mask, status, new_indices_mask);
        gpuErrchk(cudaGetLastError());

        get_num_labels(count_gpu, indices_mask, sum_scratch, M);
        gpuErrchk(cudaMemcpy(count_cpu, count_gpu, sizeof(int), cudaMemcpyDeviceToHost));
    }  

    free(count_cpu);
    cudaFree(indices_mask); cudaFree(new_indices_mask); cudaFree(known_status);
    cudaFree(count_gpu); cudaFree(minimum_gpu);
    cudaFree(min_scratch); cudaFree(sum_scratch);
}


int main() {
    int N = round((MAX - MIN) / spacing); 
    int M = N*N;
    int blockSize = 256;  // Use fixed block size for better occupancy

    // printf("%d %d", M, blockSize);

    double *solution_cpu, *solution_gpu;
    int *status_cpu, *status_gpu;

    solution_cpu = (double *) malloc (M*sizeof(double));
    status_cpu = (int *) malloc (M*sizeof(int));

    input_data(solution_cpu, status_cpu, N); 

    cudaMalloc((void**)&solution_gpu, M*sizeof(double));
    cudaMalloc((void**)&status_gpu, M*sizeof(int));

    cudaMemcpy(solution_gpu, solution_cpu, M*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(status_gpu, status_cpu, M*sizeof(int), cudaMemcpyHostToDevice);
    
    // Use reasonable grid size based on blockSize and total elements
    dim3 dimGrid((N*N + blockSize - 1) / blockSize);
    dim3 dimBlock(blockSize);

    marching_with_correctness_criterion(solution_gpu, status_gpu, dimGrid, dimBlock, M);

    cudaMemcpy(solution_cpu, solution_gpu, M*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(status_cpu, status_gpu, M*sizeof(int), cudaMemcpyDeviceToHost); 

    // real_display_array(solution_cpu, N, N);
    // int_display_array(status_cpu, N, N);  

    write_to_file_real_array(solution_cpu, N, N);

    free(solution_cpu); free(status_cpu);
    cudaFree(solution_gpu); cudaFree(status_gpu);

}