using JLD2
using SparseArrays
using PartitionedArrays

@load "matrix/matrix_coarse.jld2"  A b
A = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, A.nzval)

# Define the number of row and column partitions
num_row_partitions = 4
num_col_partitions = 4

# Partition the matrix
A_partitioned = PSparseMatrix(A, num_row_partitions, num_col_partitions)