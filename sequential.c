
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdbool.h>

#define spacing     0.0078125
#define MAX         2.0
#define MIN         0.0
#define Infi        100000.0

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
    FILE *result = fopen("sequential.txt", "w");
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

double speed_function(coordinates coors) {
    double x = coors.a;
    double y = coors.b;

    return 2 + 0.5 * sin(x) * cos(y);
}

// Convenience functions
 
coordinates indices_to_coordinates(int i, int j) {
    coordinates coors = {MIN + i*spacing, MIN + j*spacing};
    return coors;
}

bool keep_in_grid(indices idx) {

    int i,j;
    i = idx.i;
    j = idx.j;

    if (((i >= 0) && (j >= 0)) && ((i < round((MAX - MIN) / spacing)) && (j < round((MAX - MIN) / spacing)))) {
        return true;
    } 
    return false;
}

indices_list get_adjacents(int i, int j) {

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

indices_list get_diagonals(int i, int j) {

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

indices_list get_adjacents_with_diagonal(int i, int j, int p, int q) {

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

indices_list get_diagonals_with_adjacent(int i, int j, int p, int q) {

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

int get_num_label(int *indices_mask) {

    int count = 0;
    int i, j;
    int N = round((MAX - MIN) / spacing);

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (*(indices_mask + i*N + j) == 1) {
                count += 1;
            }
        }
    } 
    return count;
}

void indices_mask_init(int *indices_mask, int N) {
    int i, j;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            *(indices_mask + i*N + j) = 0;
        }
    }
}

void sum(int *indices_mask, int *new_indices_mask, int N) { 

    int i, j;
    for (i = 0; i < N; i ++) {
        for (j = 0; j < N; j++) {
            *(indices_mask + i*N + j) += *(new_indices_mask + i*N + j);
        }
    }
}

double get_min(double *solution, int *indices_mask, int N) { 

    int i, j;
    double minimum = Infi;

    for (i = 0; i < N; i ++) {
        for (j = 0; j < N; j++) {
            if (*(indices_mask + i*N + j) == 1) {
                if ((fabs(*(solution + i*N + j)) > 100001)) {
                    printf("%d %d \n", i, j);
                } 
                if (minimum > *(solution + i*N + j)) {
                    minimum = *(solution + i*N + j);
                }
            }
        }
    }
    return minimum;
}

// End convenience functions

bool correct(double *solution, int i, int j, double value, int N) {
    if (*(solution + i*N + j) < (value + (spacing / sqrt(2)))) {
        return true;
    }
    else {
        return false;
    }
}

void remove_known(int *indices_mask, int *status, int N) {

    int i, j;
    for (i = 0; i < N; i ++) {
        for (j = 0; j < N; j++) {
            if (*(status + i*N + j) == 2) {
                *(indices_mask + i*N + j) = 0;
            }
        }
    }
}

double difference_adj(double a, int i, int j) {
    return a + spacing / (speed_function(indices_to_coordinates(i,j)));
}

double difference_diag(double a, int i, int j) {
    return a + sqrt(2)*(spacing / (speed_function(indices_to_coordinates(i,j))));
}

double difference_adj_diag(double a, double b, int i, int j) {
    double difference = a - b;
    if ((difference >=0) && (difference <= (spacing / 
                                            (sqrt(2)*speed_function(indices_to_coordinates(i,j)))))) {
        return a + sqrt((pow(spacing, 2) / 
                         pow(speed_function(indices_to_coordinates(i, j)), 2)) - pow(difference,2));
    }
    return Infi;
}

void label(double *solution, int *status, int *indices_mask, int i, int j, double value) {
    int N = round((MAX - MIN) / spacing);

    if (*(status + i*N + j) == 0) {
        *(status + i*N + j) = 1;
        *(solution + i*N + j) = value;
        *(indices_mask + i*N + j) = 1;
    }
    else if ((*(status + i*N + j) == 1) && (value < *(solution + i*N + j))) {
        *(solution + i*N + j) = value;
    }
}

void label_neighbors(double *solution, int *status, int *indices_mask, int i, int j) {
    
    int N = round((MAX - MIN) / spacing);
    int k, l, p, q, m, n;
    double soln, a_soln, diag_soln;

    *(status + i*N + j) = 2;
    soln = *(solution + i*N + j);

    indices_list neighbors = get_adjacents(i,j); 
    for (k=0; k<neighbors.count; k++) {
        indices neighbor = *(neighbors.list + k); 
        p = neighbor.i;
        q = neighbor.j;
        a_soln = difference_adj(soln, p, q);

        indices_list more_neighbors = get_diagonals_with_adjacent(i,j,p,q);
        for (l = 0; l < more_neighbors.count; l++) {
            indices another_neighbor = *(more_neighbors.list + l);
            m = another_neighbor.i, n = another_neighbor.j;
            if (*(status + m*N + n) == 2) {
                a_soln = min(a_soln, difference_adj_diag(soln, *(solution + m*N + n), p, q));
            }
        } 
        free(more_neighbors.list);
        label(solution, status, indices_mask, p, q, a_soln);
    } 
    free(neighbors.list);

    neighbors = get_diagonals(i,j);
    for (k=0; k<neighbors.count; k++) {
        indices neighbor = *(neighbors.list + k);
        p = neighbor.i, q = neighbor.j;
        diag_soln = difference_diag(soln, p, q);

        indices_list more_neighbors = get_adjacents_with_diagonal(i,j,p,q);
        for (l = 0; l < more_neighbors.count; l++) {
            indices another_neighbor = *(more_neighbors.list + l); 
            m = another_neighbor.i, n = another_neighbor.j;
            if (*(status + m*N + n) == 2) {
                diag_soln = min(diag_soln, difference_adj_diag(*(solution + m*N + n), soln, p, q));
            }
        } 
        free(more_neighbors.list);
        label(solution, status, indices_mask, p, q, diag_soln);
    } 
    free(neighbors.list);
}

void marching_with_correctness_criterion(double *solution, int *status) {

    int N = round((MAX - MIN) / spacing);
    int *indices_mask, *new_indices_mask;
    int i, j, k; 

    indices_mask = (int *) malloc ((N*N)*sizeof(int));
    new_indices_mask = (int *) malloc ((N*N)*sizeof(int));    

    indices_mask_init(indices_mask, N);
    indices_mask_init(new_indices_mask, N);

    for (i = 0; i < N; i ++) {
        for (j = 0; j < N; j++) {
            if (*(status + i*N + j) == 2) { 
                label_neighbors(solution, status, indices_mask, i, j);
            }
        }
    } 

    int count = get_num_label(indices_mask);

    while (count > 0) {
        double minimum = get_min(solution, indices_mask, N);
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                if (*(indices_mask + i*N + j) == 1) {
                    if (correct(solution, i, j, minimum, N) == true) {
                        label_neighbors(solution, status, new_indices_mask, i, j);
                    }
                }
            }
        }

        remove_known(indices_mask, status, N); 
        sum(indices_mask, new_indices_mask, N);
        indices_mask_init(new_indices_mask, N);

        count = get_num_label(indices_mask);
    }
}

int main() {
    int N = round((MAX - MIN) / spacing);
    double *solution;
    int *status;
    solution = (double *) malloc (N*N*sizeof(double));
    status = (int *) malloc (N*N*sizeof(int));

    input_data(solution, status, N);

    marching_with_correctness_criterion(solution, status);

    // real_display_array(solution, N, N);
    // int_display_array(status, N, N);

    write_to_file_real_array(solution, N, N);
}