include("./Configuration.jl")

include("./FunMesh.jl")

names = ["coarse", "medium", "fine", "extrafine"]
N = [10, 15, 20, 25]

for (i, data) in enumerate(N)
    mesh_generator(names[i],
                   f, ω, L, t_P, t_F, d_PML,
                   N_rigid_circles, N_porous_circles, r_0, σ_r, tol_circle,
                   N_rigid_clamps, N_porous_clamps, r_0i, σ_ri, r_0e, σ_re, θ_o_min, θ_o_max, θ_b_min, θ_b_max, tol_clamp,
                   N_rigid_ellipses, a_0, σ_a, b_0, σ_b, e_min, e_max, θ_el_min, θ_el_max, α_min, α_max, tol_ellipse, data)
end