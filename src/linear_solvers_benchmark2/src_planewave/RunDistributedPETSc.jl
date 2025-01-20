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
      K(x) = x[2]-t_P > 0 ? K_F(ω) : K_P(ω)
      ρ(x) = x[2]-t_P > 0 ? ρ_F(ω) : ρ_P(ω)         
      γ_1(x) = abs(x[1]) < L/2  ? 1. : (x[2]-t_P < 0  ? 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2 : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2)
      γ_2(x) = x[2] > 0 ? 1. : 1. + 1im/k_P(ω) * σ_0 * x[2]^2/d_PML^2 
      H(x) = TensorValue{2,2,ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x))
      Hinv(x) = inv(H(x))
      J(x) = det(H(x))
      Jinv(x) = 1/det(H(x))
      u(x) = exact_solution_som(x, ρ_F(ω), c_F(ω), k_F(ω), ρ_P(ω), c_P(ω), k_P(ω), P_0, t_F, t_P, d_PML, σ_0)

      # Load the mesh
      model = GmshDiscreteModel(parts, "data/planewave.msh")

      # Define the finite element space: Raviart-Thomas of order 1
      order = 1
      reffe = ReferenceFE(raviart_thomas, Float64, order)
      V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides", "bottom"], vector_type=Vector{ComplexF64})

      # Define the trial function with null Dirichlet boundary conditions
      uD_null = VectorValue(0.0, 0.0)
      # uD_transducer = VectorValue(0.0, 1.0)
      U = TrialFESpace(V, [uD_null, uD_null])

      Ω = Triangulation(model) # Computational domain
      Ω_physical = Triangulation(model, tags=["physical_domain"])

      # Define the measure for the fluid and porous domains
      degree = 2
      dΩ = Measure(Ω, degree)

      # Define the measure for the transducer boundary
      Γ = BoundaryTriangulation(model, tags="transducer")
      dΓ = Measure(Γ, degree)
      n = get_normal_vector(Γ) # Normal vector to the transducer boundary

      ρ_cf = CellField(ρ, Ω)
      K_cf = CellField(K, Ω)
      H_cf = CellField(H, Ω)
      Jinv_cf = CellField(Jinv, Ω)

      # Bilinear term
      a(u, v) = ∫( K_cf * Jinv_cf * divergence(u) * divergence(v) )*dΩ -
                ∫( ρ_cf * Jinv_cf * (((H_cf) ⋅ u) ⋅ ((H_cf)  ⋅ v)) * ω^2 )*dΩ

      # Source term
      b(v) = ∫(-P_0 * (v ⋅ n))*dΓ
      

      # Assembly the system
      op = AffineFEOperator(a, b, U, V)
      # Solve the system

      # get_matrix(op)

      # options ="-pc_type gamg -ksp_type gmres"

     
      # options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart 30 -ksp_rtol 1e-5 -ksp_atol 1e-5 -ksp_converged_reason -ksp_error_if_not_converged true -ksp_monitor"
      # options = "-ksp_type cg -pc_type gamg -ksp_monitor"
      # options = "-ksp_type cg -pc_type gamg -ksp_monitor"
      # options = "-ksp_type cg -pc_type gamg -ksp_monitor"
      
      options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart 30 -ksp_atol 1.0e-10 -ksp_max_it 100000 -ksp_monitor"
      Uh = GridapPETScComplex.with(args=split(options)) do
            ls = PETScLinearSolver()
            Uh = solve(ls,op)
      end

      # solver = PETScLinearSolver()
      # Uh = solve(solver, op)
      
            # uh = (Jinv_cf) * ((Hm_cf) ⋅ Uh)
            # u = CellField(u, Ω)
            # error = u - uh

            # uh_x = (u->u[1])∘uh
            # uh_y = (u->u[2])∘uh

            # writevtk(Ω, "results/demoRT", cellfields=["Real(uh)"=>real(uh_x), "Imag(uh)"=>imag(uh_x)])

            # writevtk(Ω,"./results/result_novel"*names, cellfields=["Re(uh)"=>real(uh), "Im(uh)"=>imag(uh),
            #                                                       "Re(u)"=>real(u), "Im(u)"=>imag(u),
            #                                                       "Re(error)"=>real(error), "Im(error)"=>imag(error)])
       #GridapPetscdo
      
      # uh = (Jinv_cf) * ((Hm_cf) ⋅ Uh)
      u = CellField(u, Ω)
      # error = u - uh
      
      # uh_real = CellField(V, real(Uh.cell_dof_values))
      # uh_imag = CellField(V, imag(Uh.cell_dof_values))

      # uh_x = (u->u[1])∘uh
      uh_y = (u->u[2])∘Uh

      writevtk(Ω_physical,"./results/planewave_results/results",
               cellfields=["Re(uh)"=>real(uh_y), "Im(uh)"=>imag(uh_y),
                           "Re(uex_som)"=>real(u), "Im(u_som)"=>imag(u)])     
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