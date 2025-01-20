using PyCall
using LaTeXStrings
using PyPlot
meshio = pyimport("meshio")
rc("text",usetex="True")
rc("font",family="serif")
rc("font",size=12)

name = "translated_darcy"
labels = [L"$\mathrm{Re}(u_z)$", L"$\mathrm{Im}(u_z)$"]
files = ["Re(u_total)", "Im(u_total)"]
savers = ["./results/Re$name.svg", "./results/Im$name.svg"]
image_size_cm=(9.5, 9.5) 
size_inches = (image_size_cm[1] / 2.54, image_size_cm[2] / 2.54)

function plotVTU(name, size)
    reader = meshio.read("./results/$name.vtu")
    x = reader.points
    triangles = reader.cells_dict["triangle"]
    for i in eachindex(files)
        data = reader.point_data[files[i]]
        fig, ax = plt.subplots(figsize=size)
        contourf_plot = ax.tricontourf(x[:, 1], x[:, 2], triangles, data[:, 2], 1000, cmap="coolwarm")
        ax.set_aspect("equal")
        ax.set_rasterized(true)
        ax.set_xlabel(L"$x$")
        ax.set_ylabel(L"$z$")
        colorbar = fig.colorbar(contourf_plot, ax=ax, orientation="horizontal", pad=0.26)
        colorbar.set_label(labels[i], size=11)
        min_tick = round(minimum(data[:, 2]), sigdigits=3)
        max_tick = round(maximum(data[:, 2]), sigdigits=3)
        colorbar.set_ticks([min_tick, max_tick])
        colorbar.ax.xaxis.set_label_position("top")
        colorbar.ax.tick_params(labelsize=11)
        fig
        plt.savefig(savers[i], dpi=800, format="svg", transparent=true)
    end
end

plotVTU(name, size_inches)

