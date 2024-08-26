
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
    indices *list;
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

__global__ void mask_init(int *mask) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    *(mask + index) = 0;
}

__global__ void known_status_init(int *known_status, int *status) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
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
    int k, count = 0; 
    indices all[4] = {{i, j+1}, {i, j-1}, {i+1, j}, {i-1, j}};

    indices *list; 
    list = (indices *) malloc (4*sizeof(indices)); 

    for (k=0; k<4; k++) {
        if (keep_in_grid(all[k]) == true) { 
            *(list + count) = all[k];
            count += 1;
        }
    }

    // if (count < 4) {
    //     indices *tmp  = (indices *) realloc (list, count*sizeof(indices));
    // }

    indices_list out = {list, count};
    return out;
}

__device__ indices_list get_diagonals(indices inp_indices) {
    
    int i = inp_indices.i, j = inp_indices.j;
    int k, count = 0;
    indices all[4] = {{i+1, j+1}, {i-1, j-1}, {i+1, j-1}, {i-1, j+1}};

    indices *list;
    list = (indices *) malloc (4*sizeof(indices));

    for (k=0; k<4; k++) {
        if (keep_in_grid(all[k]) == true) {
            *(list + count) = all[k]; 

            count += 1;
        }
    }

    // if (count < 4) {
    //     indices *tmp  = (indices *) realloc (list, count*sizeof(indices));
    // }

    indices_list out = {list, count};
    return out;
}

__device__ indices_list get_adjacents_with_diagonal(indices source, indices diagonal) {
    
    int i = source.i, j = source.j;
    int p = diagonal.i, q = diagonal.j;
    int a = min(i,p), b = min(j,q);
    int k, count = 0;
    indices *list;
    list = (indices *) malloc (2*sizeof(indices));

    if (((a==i)&&(b==j)) || ((a==p)&&(b==q))) {
        indices all[2] = {{a, b+1}, {a+1, b}}; 

        for (k=0; k<2; k++) {
            if (keep_in_grid(all[k]) == true) {
                *(list + count) = all[k];
                count += 1;
            }
        }
    }
    else {
        indices all[2] = {{a, b}, {a+1, b+1}}; 

        for (k=0; k<2; k++) {
            if (keep_in_grid(all[k]) == true) {
                *(list + count) = all[k];
                count += 1;
            }
        }
    }

    // if (count < 2) {
    //     indices *tmp  = (indices *) realloc (list, count*sizeof(indices));
    // }

    indices_list out = {list, count};

    return out;
}

__device__ indices_list get_diagonals_with_adjacent(indices source, indices diagonal) {

    int i = source.i, j = source.j;
    int p = diagonal.i, q = diagonal.j;
    int k, count = 0;
    indices *list;
    list = (indices *) malloc (2*sizeof(indices));

    if (abs(i-p) > abs(j-q)) {
        indices all[2] = {{i, j+1}, {i, j-1}}; 

        for (k=0; k<2; k++) {
            if (keep_in_grid(all[k]) == true) {
                *(list + count) = all[k];
                count += 1;
            }
        }
    }
    else {
        indices all[2] = {{i+1, j}, {i-1, j}}; 

        for (k=0; k<2; k++) {
            if (keep_in_grid(all[k]) == true) {
                *(list + count) = all[k];
                count += 1;
            }
        }
    }

    // if (count < 2) {
    //     indices *tmp  = (indices *) realloc (list, count*sizeof(indices));
    // }

    indices_list out = {list, count};

    return out;
}

__global__ void sum(int *indices_mask, int *new_indices_mask) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    *(indices_mask + index) += *(new_indices_mask + index);
}

__global__ void reduce_min_1(double *solution, int *indices_mask, double *out_data) {
    int N = round((MAX - MIN) / spacing); 
    extern __shared__ double dynamic_shared_data_1[];
    // double *dynamic_shared_data = (double *)shared_data_1;
    
    int thread_idx, index;
    unsigned int s;
    thread_idx = threadIdx.x;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    if (*(indices_mask + index) == 1) {
        dynamic_shared_data_1[thread_idx] = *(solution + index);
    }
    else {
        dynamic_shared_data_1[thread_idx] = Infi;
    }
    __syncthreads();
 
    for (s = blockDim.x/2; s > 0; s>>=1) {
        if (thread_idx < s) {
            dynamic_shared_data_1[thread_idx] = min(dynamic_shared_data_1[thread_idx], dynamic_shared_data_1[thread_idx + s]);
        }
    __syncthreads();
    }

    if (thread_idx == 0) {
        *(out_data + blockIdx.x) = dynamic_shared_data_1[0];
    }
}   

__global__ void reduce_min_2(double *in_data, double *out_data) {
    extern __shared__ double dynamic_shared_data[];
    // double *dynamic_shared_data = (double *)shared_data_2;
    
    int thread_idx, index;
    unsigned int s;
    thread_idx = threadIdx.x;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    dynamic_shared_data[thread_idx] = *(in_data + index);
    __syncthreads();
 
    for (s = blockDim.x/2; s > 0; s>>=1) {
        if (thread_idx < s) {
            dynamic_shared_data[thread_idx] = min(dynamic_shared_data[thread_idx], dynamic_shared_data[thread_idx + s]);
        }
    __syncthreads();
    }

    if (thread_idx == 0) {
        *(out_data + blockIdx.x) = dynamic_shared_data[0];
    }
}

void get_min(double *minimum, double *solution, int *indices_mask, dim3 dimGrid, dim3 dimBlock) {
    int N = round((MAX - MIN) / spacing);
    double *out_data;
    gpuErrchk(cudaMalloc((void**)&out_data, GridSize*sizeof(double)));  

    size_t shareMemSize = ((N*N)/GridSize) * sizeof(double);
    reduce_min_1<<<dimGrid, dimBlock, shareMemSize>>>(solution, indices_mask, out_data);;
    gpuErrchk(cudaGetLastError());
    gpuErrchk(cudaDeviceSynchronize());
    shareMemSize = GridSize * sizeof(double);
    reduce_min_2<<<1, dimGrid, shareMemSize>>>(out_data, minimum);;
    gpuErrchk(cudaGetLastError());
    gpuErrchk(cudaDeviceSynchronize());  
}


__global__ void reduce_sum(int *in_data, int *out_data) {
    extern __shared__ int shared_data[];
    int *dynamic_shared_data = (int *)shared_data;
    
    int thread_idx, index, s;
    thread_idx = threadIdx.x;
    index = blockIdx.x * blockDim.x + threadIdx.x;
    dynamic_shared_data[thread_idx] = *(in_data + index);
    __syncthreads();
 
    for (s = blockDim.x/2; s > 0; s>>=1) {
        if (thread_idx < s) {
            dynamic_shared_data[thread_idx] += dynamic_shared_data[thread_idx + s];
        }
    __syncthreads();
    }

    if (thread_idx == 0) {
        *(out_data + blockIdx.x) = dynamic_shared_data[0];
    }
}

void get_num_labels(int *num_labels, int *indices_mask, dim3 dimGrid, dim3 dimBlock) {
    int N = round((MAX - MIN) / spacing); 
    int *out_data;
    gpuErrchk(cudaMalloc((void**)&out_data, GridSize*sizeof(int))); 

    size_t shareMemSize = ((N*N)/GridSize) * sizeof(double);
    reduce_sum<<<dimGrid, dimBlock, shareMemSize>>>(indices_mask, out_data);
    gpuErrchk(cudaGetLastError());
    gpuErrchk(cudaDeviceSynchronize());
    shareMemSize = GridSize * sizeof(double);
    reduce_sum<<<1, dimGrid, shareMemSize>>>(out_data, num_labels);;
    gpuErrchk(cudaGetLastError());
    gpuErrchk(cudaDeviceSynchronize()); 
}

__global__ void known_status_init_and_remove_known_and_sum_and_new_indices_mask_init(int *known_status, int *indices_mask, int *status, int *new_indices_mask) {
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
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
    if (correct(*(solution + index), *(minimum))) {
        *(status + index) = 2;
        *(known_status + index) = 1;
    }
}

__global__ void remove_known(int *indices_mask, int *status) {
    
    int index;
    index = blockIdx.x * blockDim.x + threadIdx.x;
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
                free(more_neighbors.list); 
                label(solution, status, indices_mask, index, a_soln);
            }
        } 
        free(neighbors.list);

        neighbors = get_diagonals(self_indices);
        for (k=0; k<neighbors.count; k++) {
            indices neighbor = *(neighbors.list + k);
            index_diag = indices_to_index(neighbor, N);

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
                free(more_neighbors.list); 
                label(solution, status, indices_mask, index, diag_soln);
            }
        } 
        free(neighbors.list);
    }
}

void marching_with_correctness_criterion(double* solution, int *status, dim3 dimGrid, dim3 dimBlock, int M) {
    int N = round((MAX - MIN) / spacing);

    int *indices_mask, *new_indices_mask, *known_status;

    gpuErrchk(cudaMalloc((void**)&indices_mask, M*sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&new_indices_mask, M*sizeof(int)));
    gpuErrchk(cudaMalloc((void**)&known_status, M*sizeof(int))); 
  
    // mask_init<<<dimGrid, dimBlock>>>(indices_mask);
    // gpuErrchk(cudaGetLastError());
    // gpuErrchk(cudaDeviceSynchronize());
    // mask_init<<<dimGrid, dimBlock>>>(new_indices_mask);
    // gpuErrchk(cudaGetLastError());
    // gpuErrchk(cudaDeviceSynchronize());
    // known_status_init<<<dimGrid, dimBlock>>>(known_status, status);
    // gpuErrchk(cudaGetLastError());
    // gpuErrchk(cudaDeviceSynchronize());

    masks_and_known_status_init<<<dimGrid, dimBlock>>>(indices_mask, new_indices_mask, known_status, status);
    gpuErrchk(cudaGetLastError());
    gpuErrchk(cudaDeviceSynchronize());

    self_label<<<dimGrid, dimBlock>>>(solution, status, known_status, indices_mask);
    gpuErrchk(cudaGetLastError());
    gpuErrchk(cudaDeviceSynchronize()); 

    mask_init<<<dimGrid, dimBlock>>>(known_status);
    gpuErrchk(cudaGetLastError());
    gpuErrchk(cudaDeviceSynchronize());

    int *count_gpu, *count_cpu;
    
    count_cpu = (int *) malloc (sizeof(int));
    gpuErrchk(cudaMalloc((void**)&count_gpu, sizeof(int)));

    get_num_labels(count_gpu, indices_mask, dimGrid, dimBlock);
    gpuErrchk(cudaMemcpy(count_cpu, count_gpu, sizeof(int), cudaMemcpyDeviceToHost));

    double *minimum_gpu;
    gpuErrchk(cudaMalloc((void**)&minimum_gpu, sizeof(double)));
 
    while (*(count_cpu) > 0) {  
        get_min(minimum_gpu, solution, indices_mask, dimGrid, dimBlock); 
        update_status_with_criterion<<<dimGrid, dimBlock>>>(solution, status, known_status, minimum_gpu);
        gpuErrchk(cudaGetLastError());
        gpuErrchk(cudaDeviceSynchronize());

        self_label<<<dimGrid, dimBlock>>>(solution, status, known_status, new_indices_mask);
        gpuErrchk(cudaGetLastError());
        gpuErrchk(cudaDeviceSynchronize());

        // mask_init<<<dimGrid, dimBlock>>>(known_status);
        // gpuErrchk(cudaGetLastError());
        // gpuErrchk(cudaDeviceSynchronize());

        // remove_known<<<dimGrid, dimBlock>>>(indices_mask, status);
        // gpuErrchk(cudaGetLastError());
        // gpuErrchk(cudaDeviceSynchronize());
        // sum<<<dimGrid, dimBlock>>>(indices_mask, new_indices_mask);
        // gpuErrchk(cudaGetLastError());
        // gpuErrchk(cudaDeviceSynchronize());
        // mask_init<<<dimGrid, dimBlock>>>(new_indices_mask);
        // gpuErrchk(cudaGetLastError());
        // gpuErrchk(cudaDeviceSynchronize());

        known_status_init_and_remove_known_and_sum_and_new_indices_mask_init<<<dimGrid, dimBlock>>>(known_status, indices_mask, status, new_indices_mask);
        gpuErrchk(cudaGetLastError());
        gpuErrchk(cudaDeviceSynchronize());

        get_num_labels(count_gpu, indices_mask, dimGrid, dimBlock);
        gpuErrchk(cudaMemcpy(count_cpu, count_gpu, sizeof(int), cudaMemcpyDeviceToHost));
    }  

    free(count_cpu);
    cudaFree(indices_mask); cudaFree(new_indices_mask); cudaFree(known_status); cudaFree(count_gpu); cudaFree(minimum_gpu);
}

int main() {
    int N = round((MAX - MIN) / spacing); 
    int M = N*N;
    int blockSize = M/GridSize;

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
    
    dim3 dimGrid(GridSize);
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


