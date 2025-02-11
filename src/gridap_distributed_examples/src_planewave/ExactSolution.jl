using LinearAlgebra


function solve_coefficients_analytical_som(ω, Z_F, Z_P, P_0, t_P, t_F, k_F, k_P)

    M = [0 0 1 0
         exp(1im*(t_P+t_F)*k_F) -exp(-1im*(t_P+t_F)*k_F) 0 0;
         exp(1im*t_P*k_F) exp(-1im*t_P*k_F) -exp(1im*t_P*k_P) -exp(-1im*t_P*k_P);
         exp(1im*t_P*k_F) -exp(-1im*t_P*k_F) -Z_P/Z_F*exp(1im*t_P*k_P) Z_P/Z_F*exp(-1im*t_P*k_P)]

    
    b = [0;
        -P_0 / (1im*ω*Z_F);
         0;
         0]

    # Resolver el sistema lineal
    solution = M \ b

    # Extraer las soluciones para A_F, B_F, A_P, y B_P
    A_F = solution[1]
    B_F = solution[2]
    A_P = solution[3]
    B_P = solution[4]

    return A_F, B_F, A_P, B_P
end

function som_check(ρ_F, c_F, k_F, ρ_P, c_P, k_P, P_0, t_F, t_P)

    # Characteristic impedance at the fluid and porous domains
    Z_F = ρ_F * c_F
    Z_P = ρ_P * c_P

    # Source term
    α = P_0/(1im*ω*Z_F)
    
    # Fluid amplitudes
    A_F = (Z_F-Z_P)/Z_F* α * exp(1im*t_F*k_F)/((Z_F+Z_P)/Z_F*exp(1im*t_P*k_F)+(Z_P-Z_F)/Z_F*exp(1im*(2*t_F+t_P)*k_F))
    B_F = A_F*exp(2im*(t_F+t_P)*k_F) + α*exp(1im*(t_F+t_P)*k_F)
    
    # Porous amplitudes
    B_P = (A_F*(exp(1im*t_P*k_F)+exp(1im*(2*t_F+t_P)*k_F))+α*exp(1im*t_F*k_F))*exp(1im*t_P*k_P) 
   
    return A_F, B_F, B_P
end

function exact_solution_som(x, ρ_F, c_F, k_F, ρ_P, c_P, k_P, P_0, t_F, t_P, d_PML, σ_0)

    # Characteristic impedance at the fluid and porous domains
    Z_F = ρ_F * c_F
    Z_P = ρ_P * c_P

    # Source term
    α = P_0/(1im*ω*Z_F)
    
    # Fluid amplitudes
    A_F = (Z_F-Z_P)/Z_F* α * exp(1im*t_F*k_F)/((Z_F+Z_P)/Z_F*exp(1im*t_P*k_F)+(Z_P-Z_F)/Z_F*exp(1im*(2*t_F+t_P)*k_F))
    B_F = A_F*exp(2im*(t_F+t_P)*k_F) + α*exp(1im*(t_F+t_P)*k_F)
    
    # Porous amplitudes
    B_P = (A_F*(exp(1im*t_P*k_F)+exp(1im*(2*t_F+t_P)*k_F))+α*exp(1im*t_F*k_F))*exp(1im*t_P*k_P) 
   
    # Evaluation at fluid, porous and bottom PML
    if (x[2] - t_P >= 0)
        return A_F * exp(1im * k_F * x[2]) + B_F * exp(-1im * k_F * x[2]) # fluid domain
    elseif (x[2] - t_P < 0 && x[2]>0)
        return B_P * exp(-1im * k_P * x[2]) # porous domain
    else
        return B_P * exp(-1im * k_P * x[2]) * exp(σ_0/3*(x[2]^3)/d_PML^2) # bottom PML (x[2]<0)
    end
end

function solve_coefficients_pml(ω, Z_F, Z_P, P_0, t_P, t_F, k_F, k_P, σ_0, d_PML)
    
    M = Complex{Float64}[
        exp(1im*(t_P+t_F)*k_F) -exp(-1im*(t_P+t_F)*k_F) 0 0 0 0;
        exp(1im*t_P*k_F) exp(-1im*t_P*k_F) -exp(1im*t_P*k_P) -exp(-1im*t_P*k_P) 0 0;
        exp(1im*t_P*k_F) -exp(-1im*t_P*k_F) -Z_P/Z_F*exp(1im*t_P*k_P) Z_P/Z_F*exp(-1im*t_P*k_P) 0 0;
        0 0 1 1 -1 -1; 
        0 0 1 -1 -1 1;  
        0 0 0 0 exp(-1im * k_P * d_PML)*exp(+d_PML*σ_0/3) exp(1im * k_P * d_PML)*exp(-d_PML*σ_0/3)]

    # Vector de términos constantes
    b = Complex{Float64}[
        -P_0 / (1im*ω*Z_F),
        0,
        0,
        0,  
        0,  
        0   
    ]

    
    solution = M \ b

    A_F = solution[1]
    B_F = solution[2]
    A_P = solution[3]
    B_P = solution[4]
    A_PML = solution[5]
    B_PML = solution[6]

    return A_F, B_F, A_P, B_P, A_PML, B_PML
end


function exact_quadratic_pml(x, t_P, k_F, k_P, d_PML, σ_0, A_F, B_F, A_P, B_P, A_PML, B_PML)
    if (x[2] - t_P >= 0)
        return A_F * exp(1im * k_F * x[2]) + B_F * exp(-1im * k_F * x[2])
    elseif (x[2] - t_P <0 && x[2]>0)
        return  A_P * exp(1im * k_P * x[2]) + B_P * exp(-1im * k_P * x[2])
    else
        return  A_PML * exp(1im * k_P * x[2]) * exp(σ_0/3*(x[2]^3)/d_PML^2) + B_PML * exp(-1im * k_P * x[2]) * exp(-σ_0/3*(-x[2]^3)/d_PML^2)
    end
end


# B_P * exp(-1im * k_P * x[2]) * exp(σ_0/3*(x[2]^3)/d_PML^2)