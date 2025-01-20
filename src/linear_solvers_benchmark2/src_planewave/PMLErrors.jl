using Gridap
using Gridap.Fields
using Gridap.Geometry
using JLD2

include("./Configuration.jl")
include("./AnalyticalSommerfeld.jl")


model = DiscreteModelFromFile("data/pml_mesh.json")

labels = get_face_labeling(model)
dimension = 2
tags = get_face_tag(labels, dimension)

fluid_PML_tag = get_tag_from_name(labels,"fluid_PML")
porous_PML_tag = get_tag_from_name(labels,"porous_PML")
bottom_PML_tag = get_tag_from_name(labels,"bottom_PML")
bottom_PML_lat_tag = get_tag_from_name(labels,"bottom_PML_lat") # Corners PML 
fluid_tag = get_tag_from_name(labels,"fluid_domain")
porous_tag = get_tag_from_name(labels,"porous_domain")

order = 1
V = TestFESpace(model, ReferenceFE(raviart_thomas, Float64, order),
    conformity=:HDiv, dirichlet_tags=["sides"], vector_type=Vector{ComplexF64})

uD = VectorValue(0.0, 0.0)

U = TrialFESpace(V, [uD])

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

k_F = 2 * pi * f/c_F
k_P = 2 * pi * f/c_P 

neumanntags = ["transducer"]
Γ = BoundaryTriangulation(model, tags=neumanntags)
dΓ = Measure(Γ, degree)
nb = get_normal_vector(Γ) # Normal vector to the transducer boundary
τ = CellField(tags, Ω) # Tags field
p = P_0 # Transducer pressure
exponents = LinRange(1, 2, 5) # De 10^1 a 10^3 con 10 puntos
logspace = exp10.(exponents)
errors = zeros(length(logspace))
aux = 1e5
sigma_best = 1
for (i, data) in enumerate(logspace)


    γ_1_F = 1 + 1im/k_F * (data)
    γ_1_P = 1 + 1im/k_P * (data)
    γ_2_P = 1 + 1im/k_P * (data) # +

    function ρ(tag)
        if tag == fluid_tag
            return ρ_F
        else
            return ρ_P
        end
    end

    function K(tag)
        if tag == fluid_tag
            return ρ_F*c_F^2
        else
            return ρ_P*c_P^2
        end
    end

    function pml_stiffness(tag, ∇u, ∇v)
        # return (TensorValue{2,2,ComplexF64}(\gamma_2(tag), 0+0im, 0+0im, \gamma_1) ⋅ (TensorValue{2,2,ComplexF64}(1/\gamma_1, 0+0im, 0+0im, 1/\gamma_2) ⊙ ∇u))     
        if tag == fluid_PML_tag
            return ((TensorValue{2,2,ComplexF64}(1/γ_1_F, 0+0im, 0+0im, 1+0im) ⊙ ∇u) ⋅ TensorValue{2,2, ComplexF64}(1+0im, 0+0im, 0+0im, γ_1_F) ⊙ ∇v)
        elseif tag == porous_PML_tag
            return ((TensorValue{2,2,ComplexF64}(1/γ_1_P, 0+0im, 0+0im, 1+0im) ⊙ ∇u) ⋅ (TensorValue{2,2, ComplexF64}(1+0im, 0+0im, 0+0im, γ_1_P) ⊙ ∇v))
        elseif tag == bottom_PML_tag
            return ((TensorValue{2,2,ComplexF64}(1+0im, 0+0im, 0+0im, 1/γ_2_P) ⊙ ∇u) ⋅ (TensorValue{2,2, ComplexF64}(γ_2_P, 0+0im, 0+0im, 1+0im) ⊙ ∇v))
        elseif tag == bottom_PML_lat_tag
            return ((TensorValue{2,2,ComplexF64}(1/γ_1_P, 0+0im, 0+0im, 1/γ_2_P) ⊙ ∇u) ⋅ (TensorValue{2,2,ComplexF64}(γ_2_P, 0+0im, 0+0im, γ_1_P) ⊙ ∇v))
        else
            return ((TensorValue{2,2,ComplexF64}(1+0im, 0+0im, 0+0im, 1+0im) ⊙ ∇u) ⋅ (TensorValue{2,2,ComplexF64}(1+0im, 0+0im, 0+0im, 1+0im) ⊙ ∇v))
        end
    end

    function pml_mass(tag)
        if tag == fluid_PML_tag
            return γ_1_F
        elseif tag == porous_PML_tag
            return γ_1_P
        elseif tag == bottom_PML_tag
            return γ_2_P
        elseif tag == bottom_PML_lat_tag
            return γ_1_P * γ_2_P
        else
            return 1.0 + 0.0im
        end
    end

    # Bilinear term
    a(u, v) = ∫( (ρ ∘ (τ)) * (pml_mass ∘ τ) * (u⋅v) * ω^2 )dΩ - ∫( (K ∘ (τ)) * (pml_stiffness ∘ (τ, ∇(u), ∇(v))) )dΩ
    # Linear temr
    b(v) = ∫( (v⋅nb) * p)dΓ
    # Assembly the system
    op = AffineFEOperator(a, b, U, V)
    # Solve the system
    uh = solve(op)
    # Post-Procesing
    uh_x = (u->u[1])∘uh
    uh_y = (u->u[2])∘uh
    function AnalyticalBox(tag) # Function to only quantifies the error in the Porous and Fluid Domains
        if tag == fluid_tag || tag == porous_tag
            return 1.0
        else
            return 0.0
        end    
    end
    # Generate a field with the tags
    λ = AnalyticalBox ∘ τ
    # Compute the analytical solution
    u = CellField(x->σ(x), Ω)
    # Compute the error and apply the mask to the error
    e = u - uh_y
    e = λ ⋅ e
    # Compute the L2 norm of
    aux2 = 100 * sqrt(sum(∫(abs2(λ*u-λ*uh_y))*dΩ)/sum(∫(abs2(λ*u))*dΩ))     
    errors[i] = aux2
    if aux2 < aux
        best_error = aux2
        global sigma_best = data
    end
    global aux = aux2
end
@save "results/errors_sigma.jld2" logspace errors 
