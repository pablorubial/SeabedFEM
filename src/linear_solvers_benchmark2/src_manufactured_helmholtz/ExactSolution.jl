function exact_solution(x, k_x, k_y)
    return exp(1im*(k_x*x[1] + k_y*x[2]))
end