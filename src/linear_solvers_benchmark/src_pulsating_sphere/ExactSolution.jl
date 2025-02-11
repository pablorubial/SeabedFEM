using SpecialFunctions

function exact_solution(x, k, x_0, y_0, r, P_0, ρ, ω)
  
    C = k * P_0 * 1/hankelh1(0, k*r) * 1/(ρ*ω^2)
    
    ux = C * hankelh1(1, k*sqrt((x[1]-x_0)^2 + (x[2]-y_0)^2)) * (x[1]-x_0) / sqrt((x[1]-x_0)^2 + (x[2]-y_0)^2)
    uy = C * hankelh1(1, k*sqrt((x[1]-x_0)^2 + (x[2]-y_0)^2)) * (x[2]-y_0) / sqrt((x[1]-x_0)^2 + (x[2]-y_0)^2)
    
    return VectorValue(ux, uy)
end
