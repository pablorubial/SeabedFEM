using Gridap
using Gridap.Geometry
using GridapGmsh
using JLD2

include("./Configuration.jl")
include("./AnalyticalIncident.jl")

function Run(freq_range, file)
    
    L2 = zeros(length(freq_range))
    # In the case that we want to save the results in an array tuples saving first frequency and second the L2 normm use this data
    # L2 = Vector{Tuple{typeof(f[1]), Float64}}(undef, length(f))
    model = GmshDiscreteModel("scenario_meshes/"*file*".msh")
    
    for (i, data) in enumerate(freq_range)

        # MeshGenerator(c_F, c_P, data, L, t_P+t_F, t_P, x_c, y_c, r_in, r_out, θ_b, θ_o, name) Not needed, take the finner mesh


        # Define the tags of the mesh boundary
        dimension = 2
        labels = get_face_labeling(model)
        tags = get_face_tag(labels, dimension)
        data
        # Define the finite element space: Raviart-Thomas of order 1
        order = 1
        reffe = ReferenceFE(raviart_thomas, Float64, order)
        V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides", "bottom", "objects"], vector_type=Vector{ComplexF64})

        # Define the trial function with null Dirichlet boundary conditions
        uD_null = VectorValue(0.0, 0.0)
        N = 100 # Number of points to discretize the interface domain [-]
        a_ = L/2 + d_PML
        
        π0 = compute_π0_coefs(N, L, d_PML, xᵦ, yᵦ, k_F(data))   
        u0, k = compute_u0_coefs(N, L, d_PML, xᵦ, yᵦ, ρ_F(data), k_F(data), data)
        π0[1] = P_0 * exp(-1im*k_F(data)*t_P) # Incident pressure field
        u0[1] = -1im * k_F(data) * P_0/(ρ_F(data)*data^2) * exp(-1im*k_F(data)*t_P) # Incident velocity field
        
        Βsel_F(kF, k) = kF > abs(k) ? 1im.*sqrt(kF^2 - k^2) : -sqrt(k^2 - kF^2)
        βF = Βsel_F.(k_F(data), k) # Bethas in the fluid domain
        Βsel_P(kP, k) = abs(kP) > abs(k) ? -1im.*sqrt(kP^2 - k^2) : +sqrt(k^2 - kP^2)
        βP = Βsel_P.(k_P(data), k) # Bethas in the porous domain
        πF_s, πP = compute_scat_coefs(π0, u0, real(ρ_F(data)), real(ρ_P(data)), βF, βP, data)
        u_incident(x) = u_incident_modes(x, t_P, t_F, L, a_, βP, πP, k, ρ_P(data), βF, πF_s, ρ_F(data), data)
        aux(x) = -u_incident(x)

        U = TrialFESpace(V, [uD_null, uD_null, aux])
        Ω = Triangulation(model) # Computational domain
        Ω_F = Triangulation(model, tags="fluid_domain_physic") 
        xp = get_physical_coordinate(Ω)

        # Define the measure for the fluid and porous domains
        degree = 2
        dΩ = Measure(Ω, degree)
        dΩ_F = Measure(Ω_F, degree)

        # Define the tensors H, H^-1 and the Jacobian for the fluid and porous PML (quadratic profile)
        K(x) = x[2]-t_P > 0 ? K_F(data) : K_P(data)
        ρ(x) = x[2]-t_P > 0 ? ρ_F(data) : ρ_P(data)         
        γ_1(x) = abs(x[1]) < L/2  ? 1. : (x[2]-t_P <= 0  ? 1. + 1im/k_P(data) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2 : 1. + 1im/k_F(data) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2)
        γ_2(x) = x[2] < 0 ? 1. + 1im/k_P(data) * σ_0 * x[2]^2/d_PML^2 : ((x[2]-(t_P+t_F)) > 0 ? 1. + 1im/k_F(data) * σ_0 * (x[2]-(t_P+t_F))^2/d_PML^2 : 1.)
        H(x) = TensorValue{2,2,ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x))
        Hinv(x) = inv(H(x))
        J(x) = det(H(x))
        Jinv(x) = 1/det(H(x))


        # Bilinear term
        a(u, v) = ∫( (K ∘ xp) * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ -
            ∫( (ρ ∘ xp) * (Jinv ∘ xp) * (((H ∘ xp) ⋅ u) ⋅ ((H ∘ xp) ⋅ v)) * data^2 )*dΩ

        # Source terms
        b(v) = 0

        # Assembly the system
        op = AffineFEOperator(a, b, U, V)

        # Solve the system
        Uh = solve(op)

        # Post-process the solution
        uh = (Jinv ∘ xp) * ((H ∘ xp) ⋅ Uh)
        u_i = u_incident ∘ xp
        ut = uh + u_i


        u_x = (u->u[1])∘ut
        u_y = (u->u[2])∘ut
    
        # In the case that we want to save the results in an array tuples saving first frequency and second the L2 normm use this data
        # L2[i] = (data, sqrt(sum(∫(abs2((u_x)))*dΩ) + sum(∫(abs2((u_y)))*dΩ)))
        L2[i] = sqrt(sum(∫(abs2((u_x)))*dΩ_F) + sum(∫(abs2((u_y)))*dΩ_F))                  
    
    end
        
    return L2     

end