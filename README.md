# A parallel approach for the Jacobi method

The final project of the Parallel and distributed systems: paradigms and model course for the A.Y. 21—22
consists in providing a parallel implementation of Jacobi method.

## Problem statement

The Jacobi iterative method computes the result of a system of equations Ax = B (with x vector of
variable of length n, A matrix of coefficients of dimension n by n and B vector of known terms of length
n) iteratively computing a new approximation of the values of the different variables according to the
formula:

$$x_i^{k+1} = \frac{1}{(A_{ii})}( bi - \sum_{j=1, j\neq{i}}^{n} A_{ij}x^k_j)$$

starting from some initial assignment of each xi (e.g. 0).
We require to implement the Jacobi method with both native C++ threads and FastFlow.

## How to compile

Run `make all` to compile the various implementations of the Jacobi method. If you want delete all files, run `make clean`

## How to run

To run executable files and launch the program, you can write in the terminal the following command:
```
./jacobi_seq n_terations dim_matrix
./jacobi_par n_iterations dim_matrix n_threads show_result
./jacobi_pinned n_iterations dim_matrix n_threads show_result
./jacobi_ff n_iterations dim_matrix n_threads show_result
```

where:
- **n_iterations**: set the number of iterations to be performed for the algorithm
- **dim_matrix**: set the size of square matrix  
- **n_threads**: set the parallelism degree 
- **show_result**: a flag to show the result of the last iteration of the algorithm. [0] no result [1] shows the result of the algorithm. 

The code will print on screen the execution time of the serial algorithm, of the parallel algorithm with both 1 and the given number of workers. Then, all the metrics computed such as speedup, efficiency and scalability are printed.
