"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using  GridapGmsh
using JLD2
using Revise

# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
includet("./AnalyticalIncident.jl")
using .AnalyticalIncident

function GenerateMatrix(name)
    
    # Load the mesh
    model = GmshDiscreteModel("data/" * "$name" * ".msh")

    # Define the finite element space: Raviart-Thomas of order 1
    order = 1
    reffe = ReferenceFE(raviart_thomas, Float64, order)
    V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides", "bottom", "objects"], vector_type=Vector{ComplexF64})

    # Define the trial function with null Dirichlet boundary conditions
    uD_null = VectorValue(0.0, 0.0)
    N = 200 # Number of points to discretize the interface domain [-]
    a_ = L/2 + d_PML
    
    π0 = compute_π0_coefs(N, L, d_PML, t_P, xᵦ, yᵦ, k_F(ω))   
    u0, k = compute_u0_coefs(N, L, d_PML, t_P, xᵦ, yᵦ, ρ_F(ω), k_F(ω), ω)
    
    Βsel_F(kF, k) = kF > abs(k) ? 1im.*sqrt(kF^2 - k^2) : -sqrt(k^2 - kF^2)
    βF = Βsel_F.(k_F(ω), k) # Bethas in the fluid domain
    Βsel_P(kP, k) = abs(kP) > abs(k) ? -1im.*sqrt(kP^2 - k^2) : +sqrt(k^2 - kP^2)
    βP = Βsel_P.(k_P(ω), k) # Bethas in the porous domain
    πF_s, πP = compute_scat_coefs(π0, u0, real(ρ_F(ω)), real(ρ_P(ω)), βF, βP, ω)
    u_incident(x) = u_incident_modes(x, t_P, t_F, L, a_, βP, πP, k, ρ_P(ω), βF, πF_s, ρ_F(ω), ω)
    aux(x) = -u_incident(x)

    U = TrialFESpace(V, [uD_null, uD_null, aux])
    Ω = Triangulation(model) # Computational domain
    xp = get_physical_coordinate(Ω)

    # Define the measure for the fluid and porous domains
    degree = 2
    dΩ = Measure(Ω, degree)

    # Define the tensors H, H^-1 and the Jacobian for the fluid and porous PML (quadratic profile)
    K(x) = x[2]-t_P > 0 ? K_F(ω) : K_P(ω)
    ρ(x) = x[2]-t_P > 0 ? ρ_F(ω) : ρ_P(ω)         
    γ_1(x) = abs(x[1]) < L/2  ? 1. : (x[2]-t_P <= 0  ? 1. + 1im/k_P(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2 : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2)
    γ_2(x) = x[2] < 0 ? 1. + 1im/k_P(ω) * σ_0 * x[2]^2/d_PML^2 : ((x[2]-(t_P+t_F)) > 0 ? 1. + 1im/k_F(ω) * σ_0 * (x[2]-(t_P+t_F))^2/d_PML^2 : 1.)
    H(x) = TensorValue{2,2,ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x))
    Hinv(x) = inv(H(x))
    J(x) = det(H(x))
    Jinv(x) = 1/det(H(x))


    # Bilinear term
    a(u, v) = ∫( (K ∘ xp) * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ -
        ∫( (ρ ∘ xp) * (Jinv ∘ xp) * (((H ∘ xp) ⋅ u) ⋅ ((H ∘ xp) ⋅ v)) * ω^2 )*dΩ

    # Source terms
    b(v) = 0

    # Assembly the system
    op = AffineFEOperator(a, b, U, V)

    # Getting the matrix and the right hand-side vector
    matrix = get_matrix(op)
    vector = get_vector(op)

    jld_filename = "matrix/matrix_" * name * ".jld2"

    matrix_name = "A"
    vector_name = "b"

    JLD2.jldopen(jld_filename, "w") do file
        file[matrix_name] = matrix  # Save the matrix with dynamic name
        file[vector_name] = vector # Save the vector with dynamic name
    end

end 

Meshes = ["coarse", "medium", "fine", "extrafine"]

for (i, data) in enumerate(Meshes)
    GenerateMatrix(data)
end

