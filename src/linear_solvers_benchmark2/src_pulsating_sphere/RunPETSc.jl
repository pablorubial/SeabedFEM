"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using  GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETScComplex
# using GridapPETSc

# Define all the neccesary functions to compute the solution
# Define the tensors H, H^-1 and the Jacobian for the fluid PML only in the horizontal direction (quadratic profile)
include("./Configuration.jl")
include("./ExactSolution.jl")

function main(distribute,nparts)
      parts = distribute(LinearIndices((prod(nparts),)))
      ω = 2π*10e3 # Frequency of the source
      # Functions to compute the tensors H, H^-1< for the fluid PML
      γ_1_aux(x, L, d_PML, σ_0, ω) = abs(x[1]) < L/2 ? 1. : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2
      γ_1(x) = γ_1_aux(x, L, d_PML, σ_0, ω)
      γ_2_aux(x, H, d_PML, σ_0, ω) = x[2] < 0 ? 1. + 1im/k_F(ω) * σ_0 * x[2]^2/d_PML^2 : ((x[2]-H) > 0 ? 1. + 1im/k_F(ω) * σ_0 * (x[2]-H)^2/d_PML^2 : 1.)
      γ_2(x) = γ_2_aux(x, H, d_PML, σ_0, ω)
      Hm(x) = TensorValue{2,2,ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x)) 
      Hinv(x) = inv(Hm(x))
      J(x) = det(Hm(x))
      Jinv(x) = 1/J(x)
      u(x) = exact_solution(x, k_F(ω), x_0, y_0, r, P_0, ρ_F(ω), ω)

      # Load the mesh
      model = GmshDiscreteModel(parts, "data/mesh.msh")

      # Define the finite element space: Raviart-Thomas of order 1
      order = 1
      reffe = ReferenceFE(raviart_thomas, Float64, order)
      V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides"], vector_type=Vector{ComplexF64})

      # Define the trial function with null Dirichlet boundary conditions
      uD_null = VectorValue(0.0, 0.0)

      U = TrialFESpace(V, [uD_null])

      degree = 2

      Ω = Triangulation(model) # Computational domain
      dΩ = Measure(Ω, degree)

      Γ = BoundaryTriangulation(model, tags="source")
      dΓ = Measure(Γ, degree)
      nb = get_normal_vector(Γ) # Normal vector to the source boundary

      Hm_cf = CellField(Hm, Ω)
      Jinv_cf = CellField(Jinv, Ω)

      # Bilinear term
      a(u, v) = ∫( ρ_F(ω) * c_F(ω)^2 * (Jinv_cf) * divergence(u) * divergence(v) )*dΩ-
                ∫( ρ_F(ω) * (Jinv_cf) * (((Hm_cf) ⋅ u) ⋅ ((Hm_cf) ⋅ v)) * ω^2 )*dΩ

      # Source term
      b(v) = ∫((v⋅nb) * P_0)dΓ # Look that the sign has been changed due to the normal vector points inside the boundary instead of outside it.

      # Assembly the system
      op = AffineFEOperator(a, b, U, V)
      # Solve the system

      # options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart 30 -ksp_rtol 1e-5 -ksp_atol 1e-5 -ksp_converged_reason -ksp_error_if_not_converged true -ksp_monitor"      
      # options = "-ksp_type cg -pc_type gamg -ksp_monitor"
      options = "-pc_type gamg -ksp_type gmres -ksp_monitor"     
      Uh = GridapPETScComplex.with(args=split(options)) do
            ls = PETScLinearSolver()
            Uh = solve(ls,op)
      end


end

nparts = (2,2)
with_mpi() do distribute
  main(distribute,nparts)
end

