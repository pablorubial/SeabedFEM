"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh
using GridapPETSc
using CairoMakie
using JLD2
using LaTeXStrings


# using GridapPETSc

# Define all the neccesary functions to compute the solution
# Define the tensors H, H^-1 and the Jacobian for the fluid PML only in the horizontal direction (quadratic profile)
include("./Configuration.jl")
include("./ExactSolution.jl")

function main(k_sweep, num, k_x, k_y)
      u(x) = exact_solution(x, k_x, k_y, L, t_P)

      # Load the mesh
      model = GmshDiscreteModel("data/manufactured_helmholtz_40.msh")

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
      options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart 30 -ksp_atol 1.0e-10 -ksp_rtol 1e-10 -ksp_max_it 10000 -ksp_monitor"

      # Solve the system and retrieve iteration statistics
      GridapPETSc.with(args=split(options)) do
            # Create the PETSc linear solver
            ls = PETScLinearSolver()
            fesolver = FESolver(ls)
            uh = zero(U)
            uh, cache = solve!(uh, fesolver, op)

            # Initialize a reference container for iteration counts
            i_petsc = Ref{PetscInt}()
            r_norm = Ref{PetscReal}()


           # Retrieve the total number of KSP iterations
            @check_error_code GridapPETSc.PETSC.KSPGetIterationNumber(cache.ksp[], i_petsc)
            @check_error_code GridapPETSc.PETSC.KSPGetResidualNorm(cache.ksp[], r_norm)
            println("Total KSP iterations: $(i_petsc[])")
            println("Residual norm: $(r_norm[])")

            e = u - uh

            # uh_real = CellField(V, real(uh.cell_dof_values))
            # uh_imag = CellField(V, imag(uh.cell_dof_values))
            u_ex = CellField(u, Ω)

            L2_err = 100 * sqrt(sum(∫(abs2(e))*dΩ))./sqrt(sum(∫(abs2(u_ex))*dΩ))
            # writevtk(Ω,"./results/results_sweep_helmholtzman/results$(num).vtu", cellfields=["Re(uh)"=>real(uh), "Im(u_h)"=>imag(uh),
            #                                                                          "Re(u)"=>real(u_ex), "Im(u)"=>imag(u_ex),
            #                                                                          "Re(error)"=>real(e), "Im(error)"=>imag(e)]) 
            return i_petsc, r_norm, L2_err
      end            
end

function _cm_to_px(cm)
      return cm * 72 / 2.54  # 72 units per inch in PDF, 2.54 cm per inch
end
  
font_path = "/usr/share/fonts/truetype/cmu/cmunrm.ttf" # Check if you have this font in your system!

# n_s = [1, 2]
# m_fr = [1, 2]
n_s = [3, 4]
m_fr = [3, 4]
colors = [:blue :green; :orange :purple]
# labels = [L"n=1, \, m=1" L"n=1, \, m=2"; L"n=2, \, m=1" L"n=2, \, m=2"]
labels = [L"n=3, \, m=3" L"n=3, \, m=4"; L"n=4, \, m=3" L"n=4, \, m=4"]

width = 16
height = 7

k = range(15, stop=30, length=100)
θs = [0, 45, 90]
# θs= 45

for θ in θs
      iterations = zeros(length(k))
      residuals = zeros(length(k))
      L2_errors = zeros(length(k))
      for (i, k_sweep) in enumerate(k)
            println("k_sweep = $k_sweep")
            k_x = k_sweep*cos(θ*pi/180)
            k_y = k_sweep*sin(θ*pi/180)
            iteration, residual, L2_error = main(k_sweep, i, k_x, k_y)
            iterations[i] = iteration[]
            residuals[i] = residual[]
            L2_errors[i] = L2_error[]
      end

      fig1 = Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
      ax1 = Axis(fig1[1, 1],
            xlabel = L"$k$",
            ylabel = "Iterations",
            # yscale = log10
            )

      lines!(ax1, k, iterations, color=:black, label= L"\theta = %$(θ)")
      # scatter!(ax1, 0,1e1, color = :red, label = "Resonances")
      for (i,n) in enumerate(n_s)
            for (j,m) in enumerate(m_fr)
                  ω_rr = c*π*sqrt((n/L)^2 + (m/t_P)^2)
                  k_rr = ω_rr/c
                  scatter!(ax1, k_rr, 1e1, color = colors[i,j], label = labels[i,j])
            end
      end
      fig1[1, 2] = Legend(fig1, ax1)
      save("./images_15_30_manufactured_fine_mesh/iterations_$(θ).pdf", fig1)

      fig2= Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
      ax2 = Axis(fig2[1, 1],
            xlabel = L"$k$",
            ylabel = L"r_k",
            # yscale = log10
            )
      lines!(ax2, k, residuals, color = :black, label = L"\theta=%$(θ)")
      # scatter!(ax2, 0, 1.0e-10, color = :red, label = "Resonances")
      for (i,n) in enumerate(n_s)
            for (j,m) in enumerate(m_fr)
                  ω_rr = c*π*sqrt((n/L)^2 + (m/t_P)^2)
                  k_rr = ω_rr/c
                  scatter!(ax2, k_rr, 1.0e-10, color = colors[i,j], label = labels[i,j])
            end
      end
      fig2[1, 2] = Legend(fig2, ax2)
      save("./images_15_30_manufactured_fine_mesh/residual_$(θ).pdf", fig2)

      fig3 = Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
      ax3 = Axis(fig3[1, 1],
            xlabel = L"$k$",
            ylabel = L"$L_2$ error",
            # yscale = log10
            )
      lines!(ax3, k, L2_errors, color = :black, label = L"\theta=%$(θ)")
      # scatter!(ax3, 0, 1e1, color = :red, label = "Resonances")
      for (i,n) in enumerate(n_s)
            for (j,m) in enumerate(m_fr)
                  ω_rr = c*π*sqrt((n/L)^2 + (m/t_P)^2)
                  k_rr = ω_rr/c
                  scatter!(ax3, k_rr, 1e1, color = colors[i,j], label = labels[i,j])
            end
      end
      fig3[1, 2] = Legend(fig3, ax3)
      save("./images_15_30_manufactured_fine_mesh/L2_error_$(θ).pdf", fig3)
end
# @save "./results/results_sweep_helmholtzman_2/results_sweep.jld2" k iterations

