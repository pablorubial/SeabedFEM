using PyPlot
rc("text",usetex="True")
rc("font",family="serif")
rc("font",size=12)
using Revise
using GridapGmsh
using Gridap
using FFTW
# Load the packages
include("Configuration.jl")
includet("AnalyticalIncidentFourier.jl")
using .AnalyticalIncidentFourier

# Construct the parameters of the problem used in the module
N = 200
n1 = 10
n2 = 12
A = [5, 5]
# Construct the parameters used in the module
params = EndFireParams(L, t_P, t_F, d_PML, Nₛ, xᵦ, yᵦ, Δy, k_F(ω), k_P(ω), σ_0, σ_0_analytical, false, ρ_F(ω), ρ_P(ω), ω, P0, A, n1, n2, N)
a_domain = params.L / 2 + params.d_PML # Half of interval length
h_fft = 2 * a_domain / N # Grid step size
x_points = range(-a_domain, stop=a_domain - h_fft, length=N) # Discretization of the interface domain [m], with size N
f1(x) = π_controlled(x, params)
# f1(x) = u_controlled(x, params)
y = f1.(x_points)

figure()
plot(x_points, real.(y), label="Real part ")
plot(x_points, imag.(y), label="Imaginary part")
axhline(0, color="black", linestyle="--")
xlabel(L"$x$ [m]")
ylabel(L"$u_z$ [m]")
legend()
tight_layout()
display(gcf())

# Compute the wavenumbers
k1 = 2 * π * fftfreq(N, 1) ./ h_fft