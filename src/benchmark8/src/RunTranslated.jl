"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl. The solution is saved in a .vtu file and
the Singled.FE.Function is saved in a .jld2 file to the compare  the results with the 
translated approach.
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using Gridap.Io
using GridapGmsh
using Revise

# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
includet("./AnalyticalIncidentFourier.jl")
using .AnalyticalIncidentFourier

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

# Define the tags of the mesh boundary
dimension = 2
labels = get_face_labeling(model)
tags = get_face_tag(labels, dimension)

# Define the finite element space: Raviart-Thomas of order 1
order = 1
reffe = ReferenceFE(raviart_thomas, Float64, order)
V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides", "bottom"], vector_type=Vector{ComplexF64})

# Define the trial function with null Dirichlet boundary conditions
uD_null = VectorValue(0.0, 0.0)
U = TrialFESpace(V, [uD_null, uD_null])

# Triangulation and spatial coordinates of the computational domain
Ω = Triangulation(model)
Ω_physical = Triangulation(model, tags=["porous_physical_domain", "fluid_physical_domain"])
ΩP = Triangulation(model, tags="porous_domain")
xp = get_physical_coordinate(Ω)
xp_physical = get_physical_coordinate(Ω_physical)


# Define the measure for the fluid and porous domains
degree = 2
dΩ = Measure(Ω, degree)
dΩP = Measure(ΩP, degree)

# Define the measure for the transducer boundary
ΓI = BoundaryTriangulation(model, tags="interface")
dΓI = Measure(ΓI, degree)
nI = get_normal_vector(ΓI) # Normal vector to the interface boundary
# Define the derived functions that arises from applied the translation of the solution technique
N = 20
n1 = 3
n2 = 4
A = [2, 2]

# Construct the parameters used in the module
params = EndFireParams(L, t_P, t_F, d_PML, Nₛ, xᵦ, yᵦ, Δy, k_F(ω), k_P(ω), σ_0, σ_0_analytical, true, ρ_F(ω), ρ_P(ω), ω, P0, A, n1, n2, N)


u_incident(x) =  uF_incident_controlled(x, params)
Π_incident(x) = πF_incident_controlled(x, params)
fₚ(x) = (ω^2 * ρ_P(ω) - k_F(ω)^2 * K_P(ω)) * u_incident(x)
g(x) = -(1-K_P(ω)/K_F(ω)) * Π_incident(x)

# Define the physical coefficients in the fluid and porous media: bulk modulus and mass density
K(x) = x[2]-t_P >= 0 ? K_F(ω) : K_P(ω)
ρ(x) = x[2]-t_P >= 0 ? ρ_F(ω) : ρ_P(ω) 

# Definition of the PML tensors H, H^-1, and the Jacobian for the fluid and porous PML (quadratic profile)
γ_1(x) = abs(x[1]) < L/2  ? 1. : (x[2]-t_P <= 0  ? 1. + 1im/k_P(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2 : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2)
γ_2(x) = x[2] < 0 ? 1. + 1im/k_P(ω) * σ_0 * x[2]^2/d_PML^2 : ((x[2]-(t_P+t_F)) > 0 ? 1. + 1im/k_F(ω) * σ_0 * (x[2]-(t_P+t_F))^2/d_PML^2 : 1.)
H(x) = TensorValue{2, 2, ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x))
Hinv(x) = inv(H(x))
J(x) = det(H(x))
Jinv(x) = 1/det(H(x))

# Bilinear term
a(u, v) = ∫( (K ∘ xp) * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ -
          ∫( (ρ ∘ xp) * (Jinv ∘ xp) * (((H ∘ xp) ⋅ u) ⋅ ((H ∘ xp) ⋅ v)) * ω^2 )*dΩ

# Source term (out always a minus symbol when the dot product with a normal vector is applied since the normal vector is pointing always inwards the surface (NumSeaHy rule))
b(v) = ∫(g * (v ⋅ nI))*dΓI + ∫(fₚ ⋅ v)*dΩP


# Assembly the FE matrices and solve the linear system
op = AffineFEOperator(a, b, U, V)
Uh = solve(op)

# Post-Processing the solution in the PML layers (uh=Uh in the physical domain)
uh = (Jinv ∘ xp) * ((H ∘ xp) ⋅ Uh)
ui = u_incident ∘ xp
ut = uh + ui

# # Save save the Singled.FE.Function in a .jld2 file
# to_jld2_file(Uh, "./results/TranslatedDarcy.jld2")

# Save the results in a vtk file and save the Singled.FE.Function in a .jld2 file
writevtk(Ω_physical,"./results/translated_darcy.vtu", cellfields=["Re(u_total)"=>real(ut), "Im(u_total)"=>imag(ut),
                                                                  "Re(u_inc)"=>real(ui), "Im(u_inc)"=>imag(ui), 
                                                                  "Re(u_scat)"=>real(uh), "Im(u_scat)"=>imag(uh)])               



