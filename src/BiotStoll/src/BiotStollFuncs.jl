module BiotStollFuncs

export compute_wave_properties, compute_wave_properties_check

using SpecialFunctions

include("SedimentModels.jl")
using ..SedimentModels: Sediment


function T(κ)
    a = exp(3 * π * 1im / 4)
    return (1/sqrt(2)*((1+1/1im)besseljx(1, κ * a)))/besseljx(0, κ * a)
end

function F(κ)
    return κ * T(κ) / (4 * (1 + 2im * T(κ) / κ))
end

function roots_bisquare_complex(a::Complex, b::Complex, c::Complex)
    # Handle trivial cases
    if a == 0 && b == 0
        if c == 0
            return "Infinite solutions"
        else
            return "No solutions"
        end
    elseif a == 0
        # Solve bx^2 + c = 0
        y = -c / b
        return [sqrt(y), -sqrt(y)]
    end
    
    # Solve the quadratic equation ay^2 + by + c = 0
    Δ = b^2 - 4 * a * c
    y1 = (-b + sqrt(Δ)) / (2 * a)
    y2 = (-b - sqrt(Δ)) / (2 * a)

    # Compute roots for x from y
    roots = []
    for y in (y1, y2)
        push!(roots, sqrt(y))
        push!(roots, -sqrt(y))
    end
    return roots
end

# Compute the wave properties of the porous domain, i.e., velocity of the longitudinal wave of first kind and attenuation coefficient
function compute_wave_properties(ω::Float64, sediment::Sediment)
    
    # Material properties with complex values
    μ = sediment.μ * (1 + 1im * sediment.δs / π)  # Shear modulus of the frame
    Kb = sediment.Kb * (1 + 1im * sediment.δl / π)  # Bulk modulus of the frame

    # Intermediate calculations
    D = sediment.Kr * (1 + sediment.β * (sediment.Kr / sediment.KF - 1))
    M = sediment.Kr^2 / (D - Kb)
    C = sediment.Kr * (sediment.Kr - Kb) / (D - Kb)
    H = (sediment.Kr - Kb)^2 / (D - Kb) + Kb + (4 / 3) * μ

    # Density and mass
    ρ = sediment.β * sediment.ρF + (1 - sediment.β) * sediment.ρr
    m = sediment.γ * sediment.ρF / sediment.β  # Added mass

    # Viscous correction term
    κ = sediment.a * sqrt(ω * sediment.ρF / sediment.η)
    Fk = F(κ)

    # Quadratic equation coefficients
    A_eq = C^2 - M * H
    B_eq = H * m * ω^2 - 1im * H * sediment.η * Fk * ω / sediment.k + ρ * ω^2 * M - 2 * sediment.ρF * ω^2 * C
    C_eq = sediment.ρF^2 * ω^4 - ρ * m * ω^4 + 1im * ρ * sediment.η * Fk * ω^3 / sediment.k

    # Solve quadratic equation for complex roots
    roots = roots_bisquare_complex(A_eq, B_eq, C_eq)
    re_root = real.(roots)
    im_root = imag.(roots)

    # Compute phase velocity (cl) and attenuation (alpha)
    cl = ω / minimum(abs.(re_root)) # [m/s]
    α = minimum(abs.(im_root)) # [Np/m]
    # In case of want to have the attenuation in dB/m multiply by 8.69
    # α = minimum(abs.(im_root)) * 8.69

    return cl, α
end

# This function is used only to check the computation of the sound velocity using directly the wavenumber in the Main.jl script
function compute_wave_properties_check(ω::Float64, sediment::Sediment)
    
    # Material properties with complex values
    μ = sediment.μ * (1 + 1im * sediment.δs / π)  # Shear modulus of the frame
    Kb = sediment.Kb * (1 + 1im * sediment.δl / π)  # Bulk modulus of the frame

    # Intermediate calculations
    D = sediment.Kr * (1 + sediment.β * (sediment.Kr / sediment.KF - 1))
    M = sediment.Kr^2 / (D - Kb)
    C = sediment.Kr * (sediment.Kr - Kb) / (D - Kb)
    H = (sediment.Kr - Kb)^2 / (D - Kb) + Kb + (4 / 3) * μ

    # Density and mass
    ρ = sediment.β * sediment.ρF + (1 - sediment.β) * sediment.ρr
    m = sediment.γ * sediment.ρF / sediment.β  # Added mass

    # Viscous correction term
    κ = sediment.a * sqrt(ω * sediment.ρF / sediment.η)
    Fk = F(κ)

    # Quadratic equation coefficients
    A_eq = C^2 - M * H
    B_eq = H * m * ω^2 - 1im * H * sediment.η * Fk * ω / sediment.k + ρ * ω^2 * M - 2 * sediment.ρF * ω^2 * C
    C_eq = sediment.ρF^2 * ω^4 - ρ * m * ω^4 + 1im * ρ * sediment.η * Fk * ω^3 / sediment.k

    # Solve quadratic equation for complex roots
    roots = roots_bisquare_complex(A_eq, B_eq, C_eq)
    
    kl_1_real = minimum(abs.(real.(roots)))
    kl_2_real = maximum(abs.(real.(roots)))
    kl_1_imag = minimum(abs.(imag.(roots)))
    kl_2_imag = maximum(abs.(imag.(roots)))

    return kl_1_real, kl_1_imag
end

### Preeliminary functions to compute the Biot-Stoll model, all thhis function needed to ocmpute the viscous correction factor F(κ) has been condensed using the scaled bessel functions that avoid the overflow problem. However, the original functions are kept here for reference

# # Kelvin function for arrays
# function kelvinb(x::AbstractArray)
#     a = exp(3 * π * 1im / 4)
#     be = besselj.(0, x .* a)  # Broadcast over x
#     ber = real.(be)           # Broadcast real part
#     bei = imag.(be)           # Broadcast imaginary part
#     return ber, bei
# end

# # Kelvin function for a single real number
# function kelvinb(x::Real)
#     a = exp(3 * π * 1im / 4)
#     be = besselj(0, x * a)
#     ber = real(be)
#     bei = imag(be)
#     return ber, bei
# end

# function kelvinb_scaled(x::Real)
#     a = exp(3 * π * 1im / 4)
#     be = besseljx(0, x * a)
#     ber = real(be)
#     bei = imag(be)
#     return ber, bei
# end

# # Kelvin derivative for arrays
# function kelvinb_derivative(x::AbstractArray)
#     a = exp(3 * π * 1im / 4)
#     be_prime = besselj.(1, x .* a)  # Broadcast over x
#     ber_prime = 1 / sqrt(2) * (real.(be_prime) + imag.(be_prime))
#     bei_prime = 1 / sqrt(2) * (imag.(be_prime) - real.(be_prime))
#     return ber_prime, bei_prime
# end

# # Kelvin derivative for a single real number
# function kelvinb_derivative(x::Real)
#     a = exp(3 * π * 1im / 4)
#     be_prime = besselj(1, x * a)
#     ber_prime = 1 / sqrt(2) * (real(be_prime) + imag(be_prime))
#     bei_prime = 1 / sqrt(2) * (imag(be_prime) - real(be_prime))
#     return ber_prime, bei_prime
# end


# function T(κ)
#     ber, bei = kelvinb(κ)
#     ber_prime, bei_prime = kelvinb_derivative(κ)  # Fixed typo here
#     return (ber_prime + 1im * bei_prime) / (ber + 1im * bei)
# end<<<<<<<<



end


