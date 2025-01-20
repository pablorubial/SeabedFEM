"""
This script allows to calculate the order of convergence of the numerical method used to solve the physical model.
"""

using CairoMakie
include("./Refine.jl")
include("./RunConvTest.jl")
h = [1, 0.5, 0.25, 0.125]
names = ["./data_convtest/coarse", "./data_convtest/medium", "./data_convtest/fine", "./data_convtest/extremelyfine"]
errors = zeros(length(names))

# Make the mesh refinement generatiing interlocked meshes
# for i in (1:(length(names)-1))
#     refine_mesh(names[i], names[i+1])
# end

# Cumpute the error for each mesh
for (i, data) in enumerate(names)
    errors[i] = Run(data*".msh")
end

# Compute the convergence order using a regression line
x = log10.(h)
y = log10.(errors)
linreg = hcat(fill!(similar(x), 1), x) \ y
m = linreg[2]
println("Convergence order is: $m")

# Save the profile of the order of convergence
f = Figure(size = (400,200), fontsize = 11)
ax = Axis(f[1, 1],
    xlabel = L"h",
    ylabel = L"100\frac{‖{u_h-u}‖_{L²(Ω)}}{‖{u}‖_{L²(Ω)}}",
    xscale = log10,
    yscale = log10,
)
lines!(ax, h, errors, color = :black)
scatter!(ax, h, errors, markersize=10, marker = :circle, color = :black)
save("./results/order_pml2.pdf", f)