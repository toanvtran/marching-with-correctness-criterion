
# Marching with Correctness Criterion for Solving Static Hamilton-Jacobi Equations in CUDA

### To run sequential code (`sequential.c`):

- Make sure your computer has `gcc`.


- In master directory:
    ```
    gcc -o runfile sequential.c -lm
    ./runfile
    ```

- Output is in `sequential.txt`

### To run parallel code (`parallel.cu`):
- Make sure your computer has `nvcc`.


- In master directory:
    ```
    nvcc -o runfile parallel.cu
    ./runfile
    ```
- Output is in `parallel.txt`

### To change input of the codes:

- Change the constants, behaviors of `input_data`, and `speed_function` in both codes.

- `N` (which equals `(MAX - MIN) / spacing`) is a positive integer.


## Reference
[Efficient Algorithms for Solving Static Hamilton-Jacobi Equations](https://thesis.library.caltech.edu/1888/)

[Optimizing Parallel Reduction in CUDA](https://imsc.uni-graz.at/haasegu/Lectures/GPU_CUDA/Lit/reduction.pdf)



