# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh
using Revise
using FFTW
# using Profile

# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
includet("./AnalyticalIncidentFourier.jl")
using .AnalyticalIncidentFourier

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

# Computational domain
Ω = Triangulation(model, tags=["porous_physical_domain", "fluid_physical_domain"]) # Computational domain
xp = get_physical_coordinate(Ω)

# Define the parameters of the problem
# Define the derived functions that arises from applied the translation of the solution technique
N = 200
n1 = 6
n2 = 10
A = [3, 5]

# Construct the parameters used in the module
params = EndFireParams(L, t_P, t_F, d_PML, Nₛ, xᵦ, yᵦ, Δy, k_F(ω), k_P(ω), σ_0, σ_0_analytical, true, ρ_F(ω), ρ_P(ω), ω, P0, A, n1, n2, N)

a = params.L / 2 + params.d_PML # Half of interval length
h_fft = 2 * a / params.N # Grid step size
# Compute the wavenumbers
k = 2 * π * fftfreq(params.N, 1) ./ h_fft
k[2]


πF_inc(x) = πF_incident_controlled(x, params)

# Interpolate that pressure field to the mesh
πF = πF_inc ∘ xp

# Write the field to a file
writevtk(Ω,"./results/controllated_presure.vtu", cellfields=["Real(πF)"=>real(πF), "Imag(πF)"=>imag(πF)])

