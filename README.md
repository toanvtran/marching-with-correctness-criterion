
# Marching with Correctness Criterion for Solving Static Hamilton-Jacobi Equations in CUDA

### To run sequential code (`sequential.c`):

- Make sure your computer have `gcc`.

- In master directory:
    ```
    gcc -o runfile sequential.c -lm
    ./runfile
    ```

- Output is in `sequential.txt`

### To run parallel code (`parallel.cu`):
- Make sure your computer have `nvcc`.

- In master directory:
    ```
    nvcc -o runfile parallel.cu
    ./runfile
    ```
- Output is in `parallel.txt`

### To change input of the codes:

- Change the constants, behaviors of `input_data`, and `speed_function` in both codes.

- Make sure that both `GridSize` and `N` (which equals `(MAX - MIN) / spacing`) in `parallel.cu` are powers of 2.


