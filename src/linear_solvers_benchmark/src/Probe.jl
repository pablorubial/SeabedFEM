using JLD2
using SparseArrays
using CairoMakie
using IterativeSolvers
using IncompleteLU

@load "matrix/matrix_coarse.jld2"  A b
A = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, A.nzval)

rows, cols, _ = findnz(A)
scatter(cols, rows, markersize=3)