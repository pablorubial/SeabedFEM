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
params = EndFireParams(L, t_P, t_F, d_PML, Nₛ, xᵦ, yᵦ, Δy, k_F(ω), k_P(ω), σ_0, σ_0_analytical, true, ρ_F(ω), ρ_P(ω), ω, P0)
N  = 200 # Number of points to discretize the interface domain [-]
a = params.L / 2 + params.d_PML # Half of interval length
h_fft = 2 * a / N # Grid step size
x_points = range(-a, stop=a - h_fft, length=N) # Discretization of the interface domain [m], with size N

# Compute the solution using the module function
f_interface(x) = π_end_fire_array(x, params)
y = f_interface.(x_points)

# Let us see how the reconstruction is performed if we try to reconstruct the function in the interface domain in a mesh with a different discretization size inherent to the finite element mesh
model = GmshDiscreteModel("data/mesh.msh")
ΓI = BoundaryTriangulation(model, tags="interface")
xp = get_physical_coordinate(ΓI)
coefs, k = compute_fft_coefs(N, y, params)

# Let us reconstruct the funciton in our interface domain by using the coordinates of the mesh
reconstruction_interface(x, coefs, k) = sum(coefs[j] * exp.(1im * k[j] * (x[1]+a)) for j in eachindex(coefs))
f_reconstruct(x) = reconstruction_interface(x, coefs, k)

# Interpolate the function on the interface domain
values = f_reconstruct ∘ xp

cell_points_interface = get_cell_points(ΓI)
points_interface = collect(cell_points_interface.cell_phys_point)
x_component = []

for i in eachindex(points_interface)
    push!(x_component, points_interface[i][1][1]) 
end
sort!(x_component)

points_interface = [Point(x, 0.4) for x in x_component]
values_on_points = evaluate(values, points_interface)

figure()
plot(x_points, real.(y), label="Real part ")
plot(x_points, imag.(y), label="Imaginary part")
plot(x_component, real.(values_on_points), label="Real part Gridap")
plot(x_component, imag.(values_on_points), label="Imaginary part Gridap")
xlabel(L"$x$ [m]")
ylabel(L"$u_z$ [m]")
legend()
tight_layout()
display(gcf())

