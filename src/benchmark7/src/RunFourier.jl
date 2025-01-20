  """
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
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
               
# Compute the solution
N = 200 # Number of points to discretize the interface domain [-]

# Construct the parameters of the problem used in the module 
params = EndFireParams(L, t_P, t_F, d_PML, Nₛ, xᵦ, yᵦ, Δy, k_F(ω), k_P(ω), σ_0, σ_0_analytical, true, ρ_F(ω), ρ_P(ω), ω, P0)

# Check anf ddisplay a warning if the maximum wavenumber reached with the FFT approximation not reaches the wavenumber of the problem
k_fft = 2 * π * fftfreq(N, 1) ./ ((L + 2*d_PML) / N) # Wavenumbers
if maximum(k_fft) < max(k_F(ω), abs(k_P(ω)))
    @warn "The maximum wavenumber reached with the FFT approximation is $(maximum(k)) and the wavenumber of the problem is $(k_F(ω))."
end 
if params.pml_active == false
    @warn "The PML is not active in the computation of the incident field, the solution may not be smooth as desirable ."
end

# Compute the scattering coefficients
πF_s, πP, βF, βP, k = compute_scat_coefs(N, params)

# Compute the scattering field (incident field in the real problem)
u_incident(x) = u_incident_modes(x, πF_s, πP, βF, βP, k, params)
u_inc(x) = u_inc_F(x, params)

# Compute the incident field (emitted field in the fluid)
# u_inc(x) = (x[2]-t_P >= 0 && abs(x[1]) <= L/2) ? VectorValue(0.0im, -1im * k_F(ω) * P0 / (ρ_F(ω) * ω^2) * exp(-1im * k_F(ω) * (x[2]-t_P))) : VectorValue(0.0 + 0.0im, 0.0 + 0.0im)

# Compute the total field

u_scat = u_incident ∘ xp 
u_inci = u_inc ∘ xp
u_total = u_scat + u_inci


# writevtk(Ω,"./results/fourier_darcy.vtu", cellfields=[ "Re(u_scat)"=>real(u_scat), "Im(u_scat)"=>imag(u_scat)])  

writevtk(Ω,"./results/fourier_darcy.vtu", cellfields=["Re(u_total)"=>real(u_total), "Im(u_total)"=>imag(u_total),
                                                      "Re(u_inc)"=>real(u_inci), "Im(u_inc)"=>imag(u_inci), 
                                                      "Re(u_scat)"=>real(u_scat), "Im(u_scat)"=>imag(u_scat)])    

                                                  

