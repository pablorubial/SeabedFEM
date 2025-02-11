function exact_solution_som(x, ρ_P, c_P, k_P, P_0, t_P, d_PML, σ_0)

    B_P = P_0/(ρ_P * c_P^2 * 1im * k_P * exp(-1im*k_P*t_P))
   
    # Evaluation at fluid, porous and bottom PML
    if ((x[2] - t_P) <= 0 && x[2] > 0)
        return B_P * exp(-1im * k_P * x[2]) # porous domain
    else
        return B_P * exp(-1im * k_P * x[2]) * exp(σ_0/3*(x[2]^3)/d_PML^2) # bottom PML (x[2]<0)
    end
end
