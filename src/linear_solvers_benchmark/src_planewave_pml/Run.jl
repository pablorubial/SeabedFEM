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
using JLD2
using CairoMakie
# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
# Include the file with the functions to compute the PML analytical solution
include("./ExactSolution.jl")

function main(ω_sweep, res, max_it)
      # Load the mesh
      model = GmshDiscreteModel("data/planewave.msh")

      # Define the tags of the mesh boundary
      dimension = 2
      labels = get_face_labeling(model)
      tags = get_face_tag(labels, dimension)

      # Define the finite element space: Raviart-Thomas of order 1
      order = 1
      reffe = ReferenceFE(raviart_thomas, Float64, order)
      V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides"], vector_type=Vector{ComplexF64})

      # Define the trial function with null Dirichlet boundary conditions
      uD_null = VectorValue(0.0, 0.0)
      # uD_transducer = VectorValue(0.0, 1.0)
      U = TrialFESpace(V, [uD_null])

      Ω = Triangulation(model) # Computational domain
      Ω_physical = Triangulation(model, tags=["physical_domain"])
      xp = get_physical_coordinate(Ω)


      # Define the measure for the fluid and porous domains
      degree = 2
      dΩ = Measure(Ω, degree)
      dΩ_physical = Measure(Ω_physical, degree)


      # Define the measure for the transducer boundary
      Γ = BoundaryTriangulation(model, tags="transducer")
      dΓ = Measure(Γ, degree)
      n = get_normal_vector(Γ) # Normal vector to the transducer boundary

      # Define the tensors H, H^-1 and the Jacobian for the fluid and porous PML (quadratic profile)        
      γ_1(x) = 1.0 + 0.0im
      γ_2(x) = x[2] > 0 ? 1.0 + 0.0im : 1. + 1im/k_P(ω_sweep) * σ_0 * x[2]^2/d_PML^2 
      H(x) = TensorValue{2,2,ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x))
      Hinv(x) = inv(H(x))
      J(x) = det(H(x))
      Jinv(x) = 1/det(H(x))
      AnalyticalBox(x) = x[2] > 0 ? 1 : 0

      # Bilinear term
      a(u, v) = ∫( K_P(ω_sweep) * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ -
            ∫( ρ_P(ω_sweep)* (Jinv ∘ xp) * (((H ∘ xp) ⋅ u) ⋅ ((H ∘ xp) ⋅ v)) * ω_sweep^2 )*dΩ

      # Source term
      b(v) = ∫(-P_0 * (v ⋅ n))*dΓ


      op = AffineFEOperator(a, b, U, V)

      # Define PETSc options adjust the restart with the number of iterations (10%)
      # Also varying maximum number of iterations (100, 1000, 10000)
      # rtol (1e-8, 1e-10, 1e-12)
      options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart $(Int(0.1*max_it)) -ksp_atol $(res) -ksp_rtol $(res) -ksp_max_it $(max_it) -ksp_monitor"

      # Solve the system and retrieve iteration statistics
      GridapPETSc.with(args=split(options)) do
            # Create the PETSc linear solver
            ls = PETScLinearSolver()
            fesolver = FESolver(ls)
            Uh = zero(U)
            Uh, cache = solve!(Uh, fesolver, op)

            # Initialize a reference container for iteration counts
            i_petsc = Ref{PetscInt}()
            r_norm = Ref{PetscReal}()

           # Retrieve the total number of KSP iterations
            @check_error_code GridapPETSc.PETSC.KSPGetIterationNumber(cache.ksp[], i_petsc)
            @check_error_code GridapPETSc.PETSC.KSPGetResidualNorm(cache.ksp[], r_norm)
            # println("Total KSP iterations: $(i_petsc[])")
            # println("Residual norm: $(r_norm[])")

            # Post process the solution
            uh = (Jinv ∘ xp) * ((H ∘ xp) ⋅ Uh)
            u_ex = CellField(x->exact_solution_som(x, ρ_P(ω_sweep), c_P(ω_sweep), k_P(ω_sweep), P_0, t_P, d_PML, σ_0), Ω)
            uh_y = (u->u[2])∘uh

            e = uh_y - u_ex

            # uh_real = CellField(V, real(uh.cell_dof_values))
            # uh_imag = CellField(V, imag(uh.cell_dof_values))

            L2_err = 100 * sqrt(sum(∫(abs2(e))*dΩ_physical))./sqrt(sum(∫(abs2(u_ex))*dΩ_physical))
            # writevtk(Ω_physical, "./results/planewave_results/results$(num).vtu", cellfields=["Re(uh)"=>real(uh_y), "Im(u_h)"=>imag(uh_y),
            #                                                                                   "Re(u)"=>real(u_ex), "Im(u)"=>imag(u_ex),
            #                                                                                   "Re(error)"=>real(e), "Im(error)"=>imag(e)]) 
            return i_petsc, r_norm, L2_err
      end            
end


function _cm_to_px(cm)
      return cm * 72 / 2.54  # 72 units per inch in PDF, 2.54 cm per inch
end
  
font_path = "/usr/share/fonts/truetype/cmu/cmunrm.ttf" # Check if you have this font in your system!
width = 16
height = 7

residual_tolerances = [1e-10, 1e-12]
max_iterations = [100, 1000, 10000]

f = range(2e3, 10e3, length=50)
ω_sweep = 2 * π * f
c_sweep = c_P.(ω_sweep)
k = ω_sweep./c_sweep

for i in residual_tolerances
      iterations = zeros(length(ω_sweep))
      residuals = zeros(length(ω_sweep))
      L2_errors = zeros(length(ω_sweep))
      for j in max_iterations
            for (k, ω) in enumerate(ω_sweep)
                  println("ω_sweep = $ω")
                  iteration, residual, L2_error = main(ω, i, j)
                  iterations[k] = iteration[]
                  residuals[k] = residual[]
                  L2_errors[k] = L2_error[]
                  @save "./results/planewave_results/results_$(j)it_res$(Int(log10(i))).jld2" iterations residuals L2_errors
            end

            fig1 = Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
            ax1 = Axis(fig1[1, 1],
                  xlabel = L"$k$",
                  ylabel = "Iterations",
                  # yscale = log10
                  )

            lines!(ax1, k, iterations, color=:black)
            save("./images_8_30_planewave_pml/iterations_$(j)it_res$(Int(log10(i))).pdf", fig1)

            fig2= Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
            ax2 = Axis(fig2[1, 1],
                  xlabel = L"$k$",
                  ylabel = L"r_k",
                  # yscale = log10
                  )
            lines!(ax2, k, residuals, color = :black)
            save("./images_8_30_planewave_pml/residual_$(j)it_res$(Int(log10(i))).pdf", fig2)

            fig3 = Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
            ax3 = Axis(fig3[1, 1],
                  xlabel = L"$k$",
                  ylabel = L"$L_2$ error",
                  # yscale = log10
                  )
            lines!(ax3, k, L2_errors, color = :black)
            save("./images_8_30_planewave_pml/L2_error_$(j)it_res$(Int(log10(i))).pdf", fig3)
      end
end

