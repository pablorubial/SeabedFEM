"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh
using GridapPETScComplex
using CairoMakie
using JLD2


# using GridapPETSc

# Define all the neccesary functions to compute the solution
# Define the tensors H, H^-1 and the Jacobian for the fluid PML only in the horizontal direction (quadratic profile)
include("./Configuration.jl")
include("./ExactSolution.jl")

function main(k_sweep, k_x, k_y)
      u(x) = exact_solution(x, k_x, k_y, L, t_P)

      # Load the mesh
      model = GmshDiscreteModel("data/manufactured_helmholtz.msh")

      # Define the finite element space: Lagrange of order 1
      order = 1
      reffe = ReferenceFE(lagrangian, Float64,order)
      V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags=["left", "right", "top", "bottom"], vector_type=Vector{ComplexF64})
      U = TrialFESpace(V, [u, u, u, u])

      Ω = Triangulation(model) # Computational domain

      # Define the measure for the fluid and porous domains
      degree = 2

      Ω = Triangulation(model)
      dΩ = Measure(Ω,degree)

      f(x) = -laplacian(x, k_x, k_y, L, t_P) - k_sweep^2*exact_solution(x, k_x, k_y, L, t_P)

      a(u,v) = ∫( (1.0+0.0im)*(∇(v)⊙∇(u)) )*dΩ - ∫( k_sweep^2*(1.0+0.0im)*(u*v) )*dΩ
      b(v) = ∫( v*f )*dΩ

      op = AffineFEOperator(a,b,U,V)

      # Define PETSc options
      options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart 30 -ksp_atol 1.0e-10 -ksp_rtol 1e-10 -ksp_max_it 100000 -ksp_monitor"

      # Solve the system and retrieve iteration statistics
      GridapPETScComplex.with(args=split(options)) do
            # Create the PETSc linear solver
            ls = PETScLinearSolver()
            fesolver = FESolver(ls)
            uh = zero(U)
            uh, cache = solve!(uh, fesolver, op)

            # Initialize a reference container for iteration counts
            i_petsc = Ref{PetscInt}()

      #     # Retrieve the total number of KSP iterations
            @check_error_code GridapPETScComplex.PETSC.KSPGetIterationNumber(cache.ksp[], i_petsc)
            println("Total KSP iterations: $(i_petsc[])")
            e = u - uh

            # uh_real = CellField(V, real(uh.cell_dof_values))
            # uh_imag = CellField(V, imag(uh.cell_dof_values))
            u_ex = CellField(u, Ω)


            # writevtk(Ω,"./results/results_sweep_helmholtzman/results$(num).vtu", cellfields=["Re(uh)"=>real(uh), "Im(u_h)"=>imag(uh),
            #                                                                          "Re(u)"=>real(u_ex), "Im(u)"=>imag(u_ex),
            #                                                                          "Re(error)"=>real(e), "Im(error)"=>imag(e)]) 
            return i_petsc
      end            
end

# k = range(0.1, stop=10, length=100)
# iterations = zeros(length(k))

# for (i, k_sweep) in enumerate(k)
#       println("k_sweep = $k_sweep")
#       θ = π/4
#       k_x = k_sweep*cos(θ)
#       k_y = k_sweep*sin(θ)
#       iteration, L2_error = main(k_sweep, k_x, k_y)
#       iterations[i] = main(k_sweep, i, k_x, k_y)[]
# end

# @save "./results/results_sweep_helmholtzman_2/results_sweep.jld2" k iterations

# function _cm_to_px(cm)
#       return cm * 72 / 2.54  # 72 units per inch in PDF, 2.54 cm per inch
# end
  
# font_path = "/usr/share/fonts/truetype/cmu/cmunrm.ttf" # Check if you have this font in your system!

# n_s = [0, 1]
# m_fr = [0, 1]

# for (i, data) in enumerate(k)
#       width = 11.2
#       height = 7
#       fig1 = Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
#       ax1 = Axis(fig1[1, 1],
#             xlabel = L"$k$",
#             ylabel = "Iterations",
#             # yscale = log10
#             )
#       lines!(ax1, k, iterations, color = :black)
#       for n in n_s
#             for m in m_fr
#                   ω_rr = c*π*sqrt((n/L)^2 + (m/t_P)^2)
#                   k_rr = ω_rr/c
#                   scatter!(ax1, k_rr, 1e1, color = :red)
#             end
#       end
#       save("./images/bench_results_100_2.pdf", fig1)
# end