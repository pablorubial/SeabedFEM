"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETScComplex
# using GridapPETSc

# Define all the neccesary functions to compute the solution
# Define the tensors H, H^-1 and the Jacobian for the fluid PML only in the horizontal direction (quadratic profile)
include("./Configuration.jl")
include("./ExactSolution.jl")

function main(distribute, nparts) #main(nparts, distribute, names, ω)
      
     
      parts = distribute(LinearIndices((prod(nparts),)))
            
      # Functions to compute the tensors H, H^-1< for the fluid PML
     # Define the tensors H, H^-1 and the Jacobian for the fluid and porous PML (quadratic profile)
      u(x) = exact_solution(x, k_x, k_y, L, t_P)

      # Load the mesh
      model = GmshDiscreteModel(parts, "data/manufactured_helmholtz.msh")

      # Define the finite element space: Raviart-Thomas of order 1
      # Define the finite element space: Lagrange of order 1
      degree = 1
      reffe = ReferenceFE(lagrangian, Float64, degree)
      V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags=["left", "right", "top", "bottom"], vector_type=Vector{ComplexF64})
      U = TrialFESpace(V, [u, u, u, u])

      Ω = Triangulation(model) # Computational domain

      # Define the measure for the fluid and porous domains
      degree = 2
      dΩ = Measure(Ω, degree)
      
      f(x) = -laplacian(x, k_x, k_y, L, t_P) - k_r^2*exact_solution(x, k_x, k_y, L, t_P)

      a(u,v) = ∫( (1.0+0.0im)*(∇(v)⊙∇(u)) )*dΩ - ∫( k_r^2*(1.0+0.0im)*(u*v) )*dΩ
      b(v) = ∫( v*f )*dΩ

      # Assembly the system
      op = AffineFEOperator(a, b, U, V)
      # Solve the system
      
      # options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart 30 -ksp_rtol 1e-5 -ksp_atol 1e-5 -ksp_converged_reason -ksp_error_if_not_converged true -ksp_monitor"
      
      options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart 30 -ksp_atol 1.0e-10 -ksp_rtol 1e-10 -ksp_max_it 100000 -ksp_monitor"
      uh = GridapPETScComplex.with(args=split(options)) do
            ls = PETScLinearSolver()
            uh = solve(ls,op)
      end
      
      u_ex = CellField(u, Ω)

      e = u_ex - uh

      writevtk(Ω,"./results/results_manufactured_bench3/results", cellfields=["Re(uh)"=>real(uh), "Im(u_h)"=>imag(uh),
                                                                  "Re(u)"=>real(u_ex), "Im(u)"=>imag(u_ex),
                                                                  "Re(error)"=>real(e), "Im(error)"=>imag(e)])     
end

# names = ["f1", "f2", "f3"]
# datas = [2*π*10e3, 2*π*12e3, 2*π*15e3]
# names = ["f1"]
# datas = [2*π*15e3]
# nparts = 2
# for i in eachindex(names)
#       with_mpi() do distribute 
#       # parts = distribute_with_mpi(LinearIndices((4,)))
#             main(nparts, distribute, names[i], datas[i])
#       end
# end

nparts = (2,2)
with_mpi() do distribute
  main(distribute, nparts)
end