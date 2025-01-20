using GridapPETScComplex: PetscScalar, PetscInt, PETSC
using SparseArrays
using JLD2

@load "matrix/matrix_coarse.jld2"  A b
# A = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, A.nzval)
# options ="-pc_type gmg -ksp_type gmres -ksp_gmres_restart 30"

# GridapPETScComplex.with(args=split(options)) do
#     solver = PETScLinearSolver()
    
  
#    end #GridapPetscdo

# A_PTSC = PETScMatrix(A)

# b_PTSC = PETScVector(b)
v_petsc = convert(PETScVector,b)
