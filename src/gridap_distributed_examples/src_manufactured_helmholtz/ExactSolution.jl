function exact_solution(x, k_x, k_y, L, t_P)
    return exp(1im*(k_x*x[1] + k_y*x[2])) * (L*x[1] - x[1]^2) * (t_P*x[2] - x[2]^2)
end

function laplacian(x, k_x, k_y, L, t_P)
    ∂²u_∂x² = (t_P*x[2]-x[2]^2) * (k_x^2*x[1]^2 - (4im*k_x + k_x^2*L)*x[1] + (2im*k_x*L - 2)) * exp(1im*(k_x*x[1] + k_y*x[2])) 
    ∂²u_∂y² = (L*x[1]-x[1]^2) * (k_y^2*x[2]^2 - (4im*k_y + k_y^2*t_P)*x[2] + (2im*k_y*t_P - 2)) * exp(1im*(k_x*x[1] + k_y*x[2]))
    return ∂²u_∂x² + ∂²u_∂y²
end