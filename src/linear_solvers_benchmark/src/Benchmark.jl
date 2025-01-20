"""
Test to use how to deal with the charge of differents matrix and delete it dynamically
"""

using JLD2
using SparseArrays
using PETSc
using LinearAlgebra


# Load data corresponding to the coarse discretisation
names = ["medium"]
size_matrix = zeros(length(names))
times_petsc = zeros(length(names))
times_LU = zeros(length(names))
error = zeros(length(names))

for (i, data) in enumerate(names)
    
    # Load the data 
    @load "matrix/matrix_" * data * ".jld2"  A b
    A = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, A.nzval)
    
    # Start the PETSc subroutine
    start_time_petsc = time()

    petsclib = PETSc.petsclibs[3]

    PETSc.initialize(petsclib)

    M = PETSc.MatSeqAIJ(A)
    println("Matrix creation completed successfully.")

    V = PETSc.VecSeq(b)
    println("Vector creation completed successfully.")

    ksp = PETSc.KSP(M; ksp_rtol=1e-8, pc_type="ilu", ksp_max_it=A.n, ksp_monitor=false)

    x_petsc = ksp \ b

    PETSc.finalize(petsclib)

    size_matrix[i] = A.n
    times_petsc[i] = time() - start_time_petsc

    # Start the LU subroutine
    start_time_LU = time()
    luA = lu(A)
    x_LU = luA\b
    times_LU[i] = time() - start_time_LU

    # Compute the error
    error[i] = norm(x_petsc - x_LU)

end




