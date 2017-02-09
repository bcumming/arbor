Aims:
* CG solver for Hines matrices arising from solution of cable equation
* written in CUDA to target GPUs
* efficiently batch process a large number matrices, which may be heterogeneous

# step 1

* Get Hines matrices
* load into test harness in C++
* write a matrix-vector multiply (CPU)

# step 2

* write CPU CG solver
* test effectiveness on Hines systems
    * sensitivity to initial guess to solution
    * efficacy of preconditioning (Jacobi and Gauss-Seidel)

# step 3

* make performance model based on these findings to understand if GPU implementation is worthwhile
    * either quite here, or continue to step 4 based on the model

# step 4

Implement basic GPU version

* basic building blocks
    * Hines gemv
    * dot product/norm
