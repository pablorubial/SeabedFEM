"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using  GridapGmsh

# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
include("./ExactSolution.jl")

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

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

U = TrialFESpace(V, [uD_null])

degree = 2

Ω = Triangulation(model) # Computational domain
dΩ = Measure(Ω, degree)

xp = get_physical_coordinate(Ω)


Γ = BoundaryTriangulation(model, tags="source")
dΓ = Measure(Γ, degree)
nb = get_normal_vector(Γ) # Normal vector to the source boundary


# Define the tensors H, H^-1 and the Jacobian for the fluid PML only in the horizontal direction (quadratic profile)
k = ω/c


# Define the tensors H, H^-1 and the Jacobian for the porous PML in the horizontal and vertical direction (quadratic profile)
γ_1(x) = abs(x[1]) < L/2 ? 1. : 1. + 1im/k * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2
γ_2(x) = x[2] < 0 ? 1. + 1im/k * σ_0 * x[2]^2/d_PML^2 : ((x[2]-H) > 0 ? 1. + 1im/k * σ_0 * (x[2]-H)^2/d_PML^2 : 1.)
Hm(x) = TensorValue{2,2,ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x)) 
Hinv(x) = inv(Hm(x))
J(x) = det(Hm(x))
Jinv(x) = 1/J(x)
AnalyticalBox(x) = (x[2] > 0 && (x[2]-H) < 0 && abs(x[1]) < L/2) ? 1 : 0

# Bilinear term
a(u, v) = ∫( ρ * c^2 * (J ∘ xp) * ((Hinv ∘ xp) ⊙ ∇(u)) ⋅ ((Hinv ∘ xp) ⊙ ∇(v)) )*dΩ - ∫( ρ * (J ∘ xp) * (u⋅v) * ω^2 )*dΩ

# a(u, v) = ∫( ρ * c^2 * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ -
      # ∫( ρ * (Jinv ∘ xp) * (((Hm ∘ xp) ⋅ u) ⋅ ((Hm ∘ xp) ⋅ v)) * ω^2 )*dΩ

# Source term
# b(v) = ∫((J ∘ xp) * P_0 *(v ⋅ nb))*dΓ
b(v) = ∫((v⋅nb) * P_0)dΓ # Look that the sign has been changed due to the normal vector points inside the boundary instead of outside it.

# Assembly the system
op = AffineFEOperator(a, b, U, V)
# Solve the system
Uh = solve(op)
uh = Uh

uh_x = (u->u[1])∘uh
uh_y = (u->u[2])∘uh

u(x) = exact_solution(x, k, x_0, y_0, r, P_0, ρ, ω)
uex = u ∘ xp
u_x = (u->u[1])∘(uex)
u_y = (u->u[2])∘(uex)

# λ = AnalyticalBox ∘ xp


writevtk(Ω,"./results/resultbad.vtu", cellfields=["Re(uh_y)"=>real(uh_y),
                                                "Im(uh_y)"=>imag(uh_y),
                                                "Im(uh_x)"=>imag(uh_x),
                                                "Re(uh_x)"=>real(uh_x),
                                                "Re(uex_y)"=>real(u_y),
                                                "Im(uex_y)"=>imag(u_y),
                                                "Im(uex_x)"=>imag(u_x),
                                                "Re(uex_x)"=>real(u_x)])

# println(100 * sqrt(sum(∫(abs2((λ*u_x-λ*uh_x)))*dΩ)/sum(∫(abs2((λ*u_x)))*dΩ)+
#                    sum(∫(abs2((λ*u_y-λ*uh_y)))*dΩ)/sum(∫(abs2((λ*u_y)))*dΩ)))