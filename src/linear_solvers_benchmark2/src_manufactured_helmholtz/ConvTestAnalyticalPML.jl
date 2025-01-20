using CairoMakie
include("./Refine.jl")
include("./RunConvTestAnalyticalPML.jl")
h = [0.5, 0.25, 0.125, 0.0625]
names = ["./data_convtest/coarse", "./data_convtest/medium", "./data_convtest/fine", "./data_convtest/extremelyfine"]
errors = zeros(length(names))

for i in (1:(length(names)-1))
    refine_mesh(names[i], names[i+1])
end

for (i, data) in enumerate(names)
    errors[i] = Run(data*".msh")
end

x = log10.(h)
y = log10.(errors)
linreg = hcat(fill!(similar(x), 1), x) \ y
m = linreg[2]
println("Convergence order is: $m")

f = Figure(size = (400,200), fontsize = 11)
ax = Axis(f[1, 1],
    xlabel = L"h",
    ylabel = L"100\frac{‖{u_h-u}‖_{L²(Ω)}}{‖{u}‖_{L²(Ω)}}",
    xscale = log10,
)
lines!(ax, h, errors, color = :black)
scatter!(ax, h, errors, markersize=10, marker = :circle, color = :black)
save("./results/order_pml_layer.pdf", f)