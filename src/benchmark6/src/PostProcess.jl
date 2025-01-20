"""
This file opens the jld2 files with the solutions of the two simulations and compares them 
computing the L2 norm of the error. It also writes the results in a .vtu file that can be
open in an unique Paraview session.
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh
using Gridap.Io
using Gridap.CellData
using Revise
# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
includet("./AnalyticalIncidentFourier.jl")
using .AnalyticalIncidentFourier

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

# Compute the solution on the fly by using the Fourier mode decomposition
N = 1 # Number of points to discretize the interface domain [-]
# Construct the parameters of the problem used in the module 
params = EndFireParams(L, t_P, t_F, d_PML, Nₛ, xᵦ, yᵦ, Δy, k_F(ω), k_P(ω), σ_0, true, ρ_F(ω), ρ_P(ω), ω, P0)
# Compute the scattering coefficients
πF_s, πP, βF, βP, k = compute_scat_coefs(N, params)
# Compute the scattering field (incident field in the real problem)
u_incident(x) = u_incident_modes(x, πF_s, πP, βF, βP, k, params)

# Define the tags of the mesh boundary
dimension = 2
labels = get_face_labeling(model)
tags = get_face_tag(labels, dimension)

# Load the Single.FE.func
Uh = from_jld2_file("./results/TranslatedDarcy.jld2") # Translated

# Get the tiangulation and the physical coordinates of the mesh
Ω = get_triangulation(model)
xp = get_physical_coordinate(Ω)
Ω_Ph = Triangulation(model, tags=["porous_physical_domain", "fluid_physical_domain"])
degree = 2
dΩ_Ph = Measure(Ω_Ph, degree) 

# Take the DOF values of the Single.FE.func == solution
vals = get_cell_dof_values(Uh)

# Fucntions to make the post-processing of the solution
γ_1(x) = abs(x[1]) < L/2  ? 1. : (x[2]-t_P <= 0  ? 1. + 1im/k_P(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2 : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2)
γ_2(x) = x[2] < 0 ? 1. + 1im/k_P(ω) * σ_0 * x[2]^2/d_PML^2 : ((x[2]-(t_P+t_F)) > 0 ? 1. + 1im/k_F(ω) * σ_0 * (x[2]-(t_P+t_F))^2/d_PML^2 : 1.)
H(x) = TensorValue{2, 2, ComplexF64}(γ_1(x), 0+0im, 0+0im, γ_2(x))
J(x) = det(H(x))
Jinv(x) = 1/det(H(x))


# Define the derived functions that arises from applied the translation of the solution technique
uF_incident(x) =  VectorValue(0.0im, -1im * k_F(ω) * P0 / (ρ_F(ω) * ω^2) * exp(-1im * k_F(ω) * (x[2]-t_P)))
u_inc(x) = (x[2]-t_P >= 0 && abs(x[1]) <= L/2) ? VectorValue(0.0im, -1im * k_F(ω) * P0 / (ρ_F(ω) * ω^2) * exp(-1im * k_F(ω) * (x[2]-t_P))) : VectorValue(0.0 + 0.0im, 0.0 + 0.0im)

# Define the finite element space previously used in perform the two simulations: Raviart-Thomas of order 1 with the corresponding boundary conditions
order = 1
reffe = ReferenceFE(raviart_thomas, Float64, order)
V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["sides", "bottom"], vector_type=Vector{ComplexF64})

# Generate the CellFields for the two solutions using the Finite Element Space and the DOF values
s1 = u_incident ∘ xp + u_inc ∘ xp # Fourier
s2 = CellField(V, vals) # Translated

# Post-Processing the solution
s2 = (Jinv ∘ xp) * ((H ∘ xp) ⋅ s2) +  uF_incident ∘ xp

# Not neccesarly to compute the L2 norm of the error, but I already had the code from another example
s1y = (u->u[2]) ∘ s1
s2y = (u->u[2]) ∘ s2

e = (s2y-s1y)/sum(∫(abs2((s2y)))*dΩ_Ph)

# Write the solution in a .vtu file
# writevtk(Ω_Ph,"./results/result_compare.vtu", cellfields=["Re(uh1)"=>real(s1), "Im(uh1)"=>imag(s1),
#                                                           "Re(uh2)"=>real(s2), "Im(uh2)"=>imag(s2),
#                                                           "Re(err)"=>real(e),  "Im(err)"=>imag(e)])     

writevtk(Ω_Ph,"./results/result_compare.vtu", cellfields=["ReFourier(uh1)"=>real(s1), "ImFourier(uh1)"=>imag(s1),
                                                          "ReTranlated(uh2)"=>real(s2), "ImTranslated(uh2)"=>imag(s2),
                                                          "Re(err)"=>real(e),  "Im(err)"=>imag(e)])

# Compute the L2 norm of the relative error
100 * sqrt(sum(∫(abs2((s2y-s1y)))*dΩ_Ph)/sum(∫(abs2((s2y)))*dΩ_Ph))                                                          