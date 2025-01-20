module AnalyticalIncidentFourier3

using Gridap
using SpecialFunctions
using FFTW


export  EndFireParams, u_end_fire_array, π_end_fire_array, u_incident_modes, compute_fft_coefs, compute_scat_coefs, u_inc_F, uF_incident, ΠF_incident

struct EndFireParams
    L::Float64
    t_P::Float64
    t_F::Float64
    d_PML::Float64
    Nₛ::Int
    xᵦ::Float64
    yᵦ::Float64
    Δy::Float64
    k_F::Float64
    k_P::ComplexF64
    σ_0::Float64
    σ_0_analytical::Float64
    pml_active::Bool
    ρ_F::Float64
    ρ_P:: ComplexF64
    ω::Float64
    P0::Float64  
end

"""
Strecthing function for the PML layer on the interface.
"""
x̂(x, params::EndFireParams) = params.pml_active ?
                             (x >= params.L/2 ?
                              x + 1im/params.k_F*params.σ_0_analytical*1/(3*params.d_PML^2) * (x - params.L/2)^3 :
                              x <= -params.L/2 ?
                              x + 1im/params.k_F*params.σ_0_analytical*1/(3*params.d_PML^2) * (x + params.L/2)^3 : x) :
                              x

"""
Derivative of the stretching function for the PML layer. 
"""
γₓ(x, params::EndFireParams) = params.pml_active ?
                               (x >= params.L/2 ?
                                1 + 1im/params.k_F*params.σ_0_analytical/ params.d_PML^2 * (x-params.L/2)^2 :
                                x <= -params.L/2 ?
                                1 + 1im/params.k_F*params.σ_0_analytical/params.d_PML^2*(x+params.L/2)^2 : 1.0 + 0.0im) :
                                1.0 + 0.0im

"""
Computation of the pressure of end-fire array on the interface between two mediums, where the Fourier mode decomposition is applied to compute the reflected and transmitted fields.
"""
function π_end_fire_array(x, params::EndFireParams)
    
    ypos = [params.yᵦ - i*params.Δy for i in 0:(params.Nₛ-1)]    
    xpos = fill(params.xᵦ, params.Nₛ)
    A = fill(1/params.Nₛ, params.Nₛ)
    
    πF = 1im/4 * sum(A[j] * hankelh1(0, params.k_F * sqrt((x̂(x,params) .- xpos[j])^2 + (params.t_P - ypos[j])^2)) for j in eachindex(ypos))
    
    return πF
end

"""
Computation of the displacement of end-fire array on the interface between two mediums, where the Fourier mode decomposition is applied to compute the reflected and transmitted fields.
"""
function u_end_fire_array(x, params::EndFireParams)
    
    ypos = [params.yᵦ - i * params.Δy for i in 0:(params.Nₛ-1)]    
    xpos = fill(params.xᵦ, params.Nₛ)
    A = fill(1/params.Nₛ, params.Nₛ)
    
    uF = -1im * params.k_F / (4 * params.ρ_F * params.ω^2) * sum(A[j] * hankelh1(1, params.k_F * sqrt((x̂(x,params) .- xpos[j])^2 + (params.t_P - ypos[j])^2)) * (params.t_P - ypos[j]) / sqrt((x̂(x,params) .- xpos[j])^2 + (params.t_P - ypos[j])^2) for j in eachindex(ypos))
    
    return uF
end

"""
Computation of the incident field in the fluid domain. This field should be added to the reflected computed field computed with Fourier Mode Decomposition to obtain the total field.
"""
function u_inc_F(x, params::EndFireParams)
    
    ypos = [params.yᵦ - i * params.Δy for i in 0:(params.Nₛ-1)]    
    xpos = fill(params.xᵦ, params.Nₛ)
    A = fill(1/params.Nₛ, params.Nₛ)
    
    if abs(x[1]) <= params.L/2 && (x[2] - params.t_P) >= 0 && (x[2] - params.t_P - params.t_F) <= 0
        uFx = -1im * params.k_F / (4 * params.ρ_F * params.ω^2) * sum(A[j] * hankelh1(1, params.k_F * sqrt((x[1] .- xpos[j])^2 + (x[2] - ypos[j])^2)) * (x[1] .- xpos[j]) / sqrt((x[1] .- xpos[j])^2 + (x[2]- ypos[j])^2) for j in eachindex(ypos))
        
        uFy = -1im * params.k_F / (4 * params.ρ_F * params.ω^2) * sum(A[j] * hankelh1(1, params.k_F * sqrt((x[1] .- xpos[j])^2 + (x[2] - ypos[j])^2)) * (x[2] - ypos[j]) / sqrt((x[1] .- xpos[j])^2 + (x[2] - ypos[j])^2) for j in eachindex(ypos))

        return VectorValue(uFx, uFy)
    else 
        return VectorValue(0.0 + 0.0im, 0.0 + 0.0im)
    end
end

function uF_incident(x, params::EndFireParams)
    
    ypos = [params.yᵦ - i * params.Δy for i in 0:(params.Nₛ-1)]    
    xpos = fill(params.xᵦ, params.Nₛ)
    A = fill(1/params.Nₛ, params.Nₛ)
    
    uFx = -1im * params.k_F / (4 * params.ρ_F * params.ω^2) * sum(A[j] * hankelh1(1, params.k_F * sqrt((x[1] .- xpos[j])^2 + (x[2] - ypos[j])^2)) * (x[1] .- xpos[j]) / sqrt((x[1] .- xpos[j])^2 + (x[2]- ypos[j])^2) for j in eachindex(ypos))
    
    uFy = -1im * params.k_F / (4 * params.ρ_F * params.ω^2) * sum(A[j] * hankelh1(1, params.k_F * sqrt((x[1] .- xpos[j])^2 + (x[2] - ypos[j])^2)) * (x[2] - ypos[j]) / sqrt((x[1] .- xpos[j])^2 + (x[2] - ypos[j])^2) for j in eachindex(ypos))

    return VectorValue(uFx, uFy)
   
end

function ΠF_incident(x, params::EndFireParams)
    
    ypos = [params.yᵦ - i * params.Δy for i in 0:(params.Nₛ-1)]    
    xpos = fill(params.xᵦ, params.Nₛ)
    A = fill(1/params.Nₛ, params.Nₛ)
    
    ΠF = 1im/4 * sum(A[j] * hankelh1(0, params.k_F * sqrt((x[1] .- xpos[j])^2 + (x[2] - ypos[j])^2)) for j in eachindex(ypos))

    return ΠF
   
end

function compute_fft_coefs(N::Int, func::Function, params::EndFireParams)
    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / N # Grid step size
    x_points = range(-a, stop=a - h_fft, length=N) # Discretization of the interface domain [m], with size N

    # Compute the wavenumbers
    k = 2 * π * fftfreq(N, 1)./h_fft

    # Compute the values of the function on the grid points
    f_values = map(x -> func(x, params), x_points)

    # Compute the Fourier coefficients of the input data and normalize by N using the FFT
    coefs = fft(f_values) ./ N

    return coefs, k, h_fft
end


function compute_fft_coefs(N::Int, points::AbstractVector, params::EndFireParams)
    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / N # Grid step size

    # Compute the wavenumbers
    k = 2 * π * fftfreq(N, h_fft)./h_fft

    # Compute the Fourier coefficients of the input data and normalize by N using the FFT
    coefs = fft(points) ./ N

    return coefs, k, h_fft
end

function compute_scat_coefs(N, params::EndFireParams)
    
    π0, k, h_fft = compute_fft_coefs(N, π_end_fire_array, params)
    u0, _ = compute_fft_coefs(N, u_end_fire_array, params)

    
    # π0[1] = params.P0 # Incident pressure field
    # u0[1] = -1im * params.k_F * params.P0/(params.ρ_F*params.ω^2) # Incident displacement field

    Βsel_F(k_F, k) = k_F > abs(k) ? 1im.*sqrt(k_F^2 - k^2) : -sqrt(k^2 - k_F^2)
    Βsel_P(k_P, k) = abs(k_P) > abs(k) ? -1im.*sqrt(k_P^2 - k^2) : +sqrt(k^2 - k_P^2)
    
    βF = Βsel_F.(params.k_F, k) # Bethas in the fluid domain
    βP = Βsel_P.(params.k_P, k) # Bethas in the porous domain
    
    # Coefficients of the 2x2 linear system associated with the coupling conditions on the interface
    A = 1 / (params.ρ_F * params.ω^2) .* βF
    B = 1 / (params.ρ_P * params.ω^2) .* βP

    # Solution of the 2x2 linear system
    πF_s = (B .* π0 .- u0) ./ (A .- B)
    πP = πF_s .+ π0

    return πF_s, πP, βF, βP, k, h_fft
end

# function uF_s(x, βF, πF_s, k, h_fft, params::EndFireParams)
    
#     ux = 1 / (params.ρ_F * params.ω^2) * sum(πF_s[j] * 1im * k[j] * exp(βF[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + params.L/2 + params.d_PML)/h_fft) for j in eachindex(βF))
    
#     uy = 1 / (params.ρ_F * params.ω^2) * sum(πF_s[j] * βF[j] * exp(βF[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + params.L/2 + params.d_PML)/h_fft) for j in eachindex(βF))
    
#     return VectorValue(ux, uy)
# end

function uF_s(x, βF, πF_s, k, h_fft, params::EndFireParams)
    
    ux = 1 / (params.ρ_F * params.ω^2) * sum(πF_s[j] * 1im * k[j] *  exp(βF[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + params.L/2 + params.d_PML)) for j in eachindex(βF))
    
    uy = 1 / (params.ρ_F * params.ω^2) * sum(πF_s[j] * βF[j] * exp(βF[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + params.L/2 + params.d_PML)) for j in eachindex(βF))
    
    return VectorValue(ux, uy)
end

# function uP(x, βP, πP, k, h_fft, params::EndFireParams)
    
#     ux = 1 / (params.ρ_P * params.ω^2) * sum(πP[j] * 1im * k[j] * exp(βP[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + params.L/2 + params.d_PML)/h_fft) for j in eachindex(βP))
    
#     uy = 1 / (params.ρ_P * params.ω^2) * sum(πP[j] * βP[j] * exp(βP[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + params.L/2 + params.d_PML)/h_fft) for j in eachindex(βP))
    
#     return VectorValue(ux, uy)
# end

function uP(x, βP, πP, k, h_fft, params::EndFireParams)
    
    ux = 1 / (params.ρ_P * params.ω^2) * sum(πP[j] * 1im * k[j] * exp(βP[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + params.L/2 + params.d_PML)) for j in eachindex(βP))
    
    uy = 1 / (params.ρ_P * params.ω^2) * sum(πP[j] * βP[j] * exp(βP[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + params.L/2 + params.d_PML)) for j in eachindex(βP))
    
    return VectorValue(ux, uy)
end

function u_incident_modes(x, πF_s, πP, βF, βP, k, h_fft, params::EndFireParams)    
    
    if abs(x[1]) <= params.L/2 && (x[2] - params.t_P) < 0 && x[2] >= 0
        return uP(x, βP, πP, k, h_fft, params) 
    elseif abs(x[1]) <= params.L/2 && (x[2] - params.t_P) >= 0 && (x[2] - (params.t_P + params.t_F)) <= 0
        return uF_s(x, βF, πF_s, k, h_fft, params)
    else
        return VectorValue(0.0 + 0.0im, 0.0 + 0.0im)
    end
end

end