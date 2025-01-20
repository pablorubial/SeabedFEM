# Script to check the real and imaginary part of the sound velocity for a sediment choosen in the Configuration.jl file
using CairoMakie
using LaTeXStrings

include("Configuration.jl")


f = [10^i for i in range(0, 6, length=10000)]
ω = 2π .* f
CP = C_P.(ω)


function _cm_to_px(cm)
    return cm * 72 / 2.54  # 72 units per inch in PDF, 2.54 cm per inch
end


# font_path = "/usr/share/texmf/fonts/opentype/public/tex-gyre-math/texgyreschola-math.otf"
font_path = "/usr/share/fonts/truetype/cmu/cmunrm.ttf"
## Plot the results in log-scale
sediments = ["MediumSilt", "VeryFineSand", "MediumSand"]

for (i, data) in enumerate(sediments)
    width = 11.2
    height = 7
    sediment = predefined_sediment(data; ρF=ρ_F(ω), KF=K_F(ω), η=η_F)
    func(f) = compute_wave_properties(f, sediment)
    freqs = 10 .^ range(0, 6, length=10000)
    res = func.(freqs)
    cl = [r[1] for r in res]
    α = [r[2] for r in res] .*8.69
    fig1 = Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
    ax1 = Axis(fig1[1, 1],
        xscale = log10,
        xlabel = L"$f$ [Hz]",
        ylabel = L"$c_l$ [m/s]",
        )
    lines!(ax1, freqs, cl, label = L"$c_l$", color = :black)
    xlims!(ax1, 1e0, 1e6)
    ylims!(ax1, 1400, 1900)
    fig2 = Figure(size = (_cm_to_px(width), _cm_to_px(height)), fontsize = 12, fonts = (; regular = font_path))
    save("./doc/figures/fig1_cl_$data.pdf", fig1)
    ax2 = Axis(fig2[1, 1],
        xscale = log10,
        yscale = log10,
        xlabel = L"$f$ [Hz]",
        ylabel = L"$\alpha$ [dB/m]",
        )
    lines!(ax2, freqs, α, label = L"$\alpha$", color = :black)
    xlims!(ax2, 1e0, 1e6)
    ylims!(ax2, 1e-5, 1e3)
    save("./doc/figures/fig2_alpha_$data.pdf", fig2)
end
