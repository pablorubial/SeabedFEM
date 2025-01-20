module AnalyticalIncidentFourier

using Gridap
using SpecialFunctions
using FFTW


export  EndFireParams, u_controlled, π_controlled, u_incident_modes, compute_fft_coefs, u_inc_F, uF_incident, ΠF_incident, πF_incident_controlled, uF_incident_controlled, u_inc_F_controlled, compute_scat_coefs


struct EndFireParams
    L::Float64 # Horizontal length of the domain [m]
    t_P::Float64 # Thickness of the porous domain [m]
    t_F::Float64 # Thickness of the fluid domain [m]
    d_PML::Float64 # Thickness of the PML layer [m]
    Nₛ::Int # Number of points to simulate the sonar transducer seen as Dirac's deltas [-]
    xᵦ::Float64 # Horizontal coordinate of the sonar transducer [m]
    yᵦ::Float64 # Vertical coordinate of the sonar transducer [m]
    Δy::Float64 # Separation between the points of the sonar discretization points [m]
    k_F::Float64 # Wavenumber in the fluid domain [rad/m]
    k_P::ComplexF64 # Wavenumber in the porous domain [rad/m]
    σ_0::Float64 # Absorption coefficient of the PML layer in the FEM simulation 
    σ_0_analytical::Float64 # Absorption coefficient of the PML layer in the analytical solution to use Fourier Mode Decomposition
    pml_active::Bool # Flag to activate the PML layer
    ρ_F::Float64 # Density of the fluid domain [kg/m³]
    ρ_P:: ComplexF64 # Density of the porous domain [kg/m³]
    ω::Float64 # Angular frequency of the sonar transducer [rad/s]
    P0::Float64 # Transducer pressure [Pa]
    A::AbstractVector # Amplitude correspondent to each of the nodes
    n1::Int # Index of the first mode chosen
    n2::Int # Index of the second mode chosen
    N::Int # Number of basis function to use in the Fourier Mode Decomposition
end

"""
Strecthing function for the PML layer on the interface.
"""
x̂(x, params::EndFireParams) = params.pml_active ?
                             (x >= params.L/2 ?
                              x + 1im/params.k_F*params.σ_0_analytical*1/(3*params.d_PML^2) * (x - params.L/2)^3 :
                              x <= -params.L/2 ?
                              x - 1im/params.k_F*params.σ_0_analytical*1/(3*params.d_PML^2) * (x + params.L/2)^3 :
                              x) :
                              x

"""
Derivative of the stretching function for the PML layer. It is not necceary to use this function in the current implementation.
"""
γₓ(x, params::EndFireParams) = params.pml_active ?
                               (x >= params.L/2 ?
                                1 + 1im/params.k_F*params.σ_0_analytical/params.d_PML^2 * (x-params.L/2)^2 :
                                x <= -params.L/2 ?
                                1 + 1im/params.k_F*params.σ_0_analytical/params.d_PML^2 * (x+params.L/2)^2 :
                                1.0 + 0.0im) :
                                1.0 + 0.0im

"""
Pressure field to be used in the Fourier Mode Decomposition. Note that the z-component is evaluated on the interface to make the Fourier Mode Decomposition.
"""
function π_controlled(x, params::EndFireParams)
   
    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / params.N # Grid step size
    # Compute the wavenumbers
    k = 2 * π * fftfreq(params.N, 1) ./ h_fft
    k₁ = k[params.n1]
    k₂ = k[params.n2]

    πF = params.A[1] * exp(1im * k₁ * x̂(x, params)) * exp(-1im * k₁ * params.t_P) + 
         params.A[2] * exp(1im * k₂ * x̂(x, params)) * exp(-1im * k₂ * params.t_P)
    
    return πF

end

"""
Incident field to be used in the Fourier Mode Decomposition. Note that the z-component is evaluated on the interface to make the Fourier Mode Decomposition.
"""
function u_controlled(x, params::EndFireParams)

    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / params.N # Grid step size
    # Compute the wavenumbers
    k = 2 * π * fftfreq(params.N, 1) ./ h_fft
    k₁ = k[params.n1]
    k₂ = k[params.n2]

    uF = 1 / (params.ρ_F * params.ω^2) * (
        - params.A[1] * k₁ * exp(1im * k₁ * x̂(x, params)) * exp(-1im * k₁ * params.t_P)
        - params.A[2] * k₂ * exp(1im * k₂ * x̂(x, params)) * exp(-1im * k₂ * params.t_P)
        )
     
    return uF

end

"""
Computation of the incident field in the fluid domain. This field should be added to the reflected computed field using Fourier Mode Decomposition to obtain the total field.
"""
function u_inc_F_controlled(x, params::EndFireParams)
    
    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / params.N # Grid step size
    # Compute the wavenumbers
    k = 2 * π * fftfreq(params.N, 1) ./ h_fft
    k₁ = k[params.n1]
    k₂ = k[params.n2]
    if abs(x[1]) <= params.L/2 && (x[2] - params.t_P) >= 0 && (x[2] - params.t_P - params.t_F) <= 0
        
        ux = 1 / (params.ρ_F * params.ω^2) * (
            params.A[1] * 1im * k₁ * exp(1im * k₁ * x[1]) * exp(-1im * k₁ * x[2]) + 
            params.A[2] * 1im * k₂ * exp(1im * k₂ * x[1]) * exp(-1im * k₂ * x[2])
            )
        
        uy = 1 / (params.ρ_F * params.ω^2) * (
            -params.A[1] * k₁ * exp(1im * k₁ * x[1]) * exp(-1im * k₁ * x[2])
            -params.A[2] * k₂ * exp(1im * k₂ * x[1]) * exp(-1im * k₂ * x[2]))

        return VectorValue(ux, uy)
    else 
        return VectorValue(0.0 + 0.0im, 0.0 + 0.0im)
    end
end

"""
Displacements incident field in all the domains. This field is used in the translation of the solution technique as incident displacement field in porous and fluid domain to compute the solution through the Finite Element Method.
"""
function uF_incident_controlled(x, params::EndFireParams)
    
    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / params.N # Grid step size
    # Compute the wavenumbers
    k = 2 * π * fftfreq(params.N, 1) ./ h_fft
    k₁ = k[params.n1]
    k₂ = k[params.n2]

    ux = 1 / (params.ρ_F * params.ω^2) * (
        params.A[1] * 1im * k₁ * exp(1im * k₁ * x[1]) * exp(-1im * k₁ * x[2]) + 
        params.A[2] * 1im * k₂ * exp(1im * k₂ * x[1]) * exp(-1im * k₂ * x[2])
        )
    
    uy = 1 / (params.ρ_F * params.ω^2) * (
        -params.A[1] * k₁ * exp(1im * k₁ * x[1]) * exp(-1im * k₁ * x[2])
        -params.A[2] * k₂ * exp(1im * k₂ * x[1]) * exp(-1im * k₂ * x[2]))

    return VectorValue(ux, uy)

end

"""
Pressure incident field in all the domains. This field is used in translation of the solution technique as incident pressure field in porous and fluid domain to compute the solution through the Finite Element Method.
"""
function πF_incident_controlled(x, params::EndFireParams)
    
    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / params.N # Grid step size
    # Compute the wavenumbers
    k = 2 * π * fftfreq(params.N, 1) ./ h_fft
    k₁ = k[params.n1]
    k₂ = k[params.n2]
    # πF = params.A[1] * exp(1im * k₁ * x[1]) * exp(-1im * 2*k₁ * x[2]) + 
    #      params.A[2] * exp(1im * k₂ * x[1]) * exp(-1im * 2*k₂ * x[2])
    πF = params.A[1] * exp(-1im * (k₁ * x[1] + 2*k₁ * x[2])) + 
         params.A[2] * exp(-1im * (k₂ * x[1] + 2*k₂ * x[2]))

    return πF

end

"""
This function computes the FFT using N basis function of the function func evaluated on the interface between the two domains. (z=t_P)
"""
function compute_fft_coefs(N::Int, func::Function, params::EndFireParams)
    a = params.L / 2 + params.d_PML # Half of interval length
    h_fft = 2 * a / N # Grid step size
    x_points = range(-a, stop=a - h_fft, length=N) # Discretization of the interface domain [m], with size N

    # Compute the wavenumbers
    k = 2 * π * fftfreq(N, 1) ./ h_fft

    # Compute the values of the function on the grid points
    f_values = map(x -> func(x, params), x_points)

    # Compute the Fourier coefficients of the input data and normalize by N using the FFT
    coefs = fft(f_values) ./ N

    return coefs, k
end


# function compute_fft_coefs(N::Int, points::AbstractVector, params::EndFireParams)
#     a = params.L / 2 + params.d_PML # Half of interval length
#     h_fft = 2 * a / N # Grid step size

#     # Compute the wavenumbers
#     k = 2 * π * fftfreq(N, 1) ./ h_fft

#     # Compute the Fourier coefficients of the input data and normalize by N using the FFT
#     coefs = fft(points) ./ N

#     return coefs, k
# end

"""
This function computes the scattering coefficients used for construct the solution using Fourier modes.
"""
function compute_scat_coefs(N, params::EndFireParams)
    
    π0, k = compute_fft_coefs(N, π_controlled, params)
    u0, _ = compute_fft_coefs(N, u_controlled, params)

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

    return πF_s, πP, βF, βP, k
end

"""
Construction of the reflected field in the fluiid domain using the scattering coefficients computed. 
"""
function uF_s(x, βF, πF_s, k, params::EndFireParams)
    
    a = params.L/2 + params.d_PML
    
    ux = 1 / (params.ρ_F * params.ω^2) * sum(πF_s[j] * 1im * k[j] * exp(βF[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + a)) for j in eachindex(βF))
    
    uy = 1 / (params.ρ_F * params.ω^2) * sum(πF_s[j] * βF[j] * exp(βF[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + a)) for j in eachindex(βF))
    
    return VectorValue(ux, uy)
end

"""
Construction of the transmitted field in the porous domain using the scattering coefficients computed.
"""
function uP(x, βP, πP, k, params::EndFireParams)

    a = params.L/2 + params.d_PML

    ux = 1 / (params.ρ_P * params.ω^2) * sum(πP[j] * 1im * k[j] * exp(βP[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + a)) for j in eachindex(βP))
    
    uy = 1 / (params.ρ_P * params.ω^2) * sum(πP[j] * βP[j] * exp(βP[j] * (x[2] - params.t_P)) * exp(1im * k[j] * (x[1] + a)) for j in eachindex(βP))
    
    return VectorValue(ux, uy)
end

"""
Combination of the two previous fields, this function is the one used in Gridap.
"""
function u_incident_modes(x, πF_s, πP, βF, βP, k,  params::EndFireParams)    
    
    if abs(x[1]) <= params.L/2 && (x[2] - params.t_P) < 0 && x[2] >= 0
        return uP(x, βP, πP, k, params) 
    elseif abs(x[1]) <= params.L/2 && (x[2] - params.t_P) >= 0 && (x[2] - (params.t_P + params.t_F)) <= 0
        return uF_s(x, βF, πF_s, k, params)
    else
        return VectorValue(0.0 + 0.0im, 0.0 + 0.0im)
    end
end

end