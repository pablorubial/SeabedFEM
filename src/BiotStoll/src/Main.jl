# Script to test if the wave properties computed splitting the computation of the wavenumbers and the introdtuction of the attenuation coefficient are correct
using LinearAlgebra
include("Configuration.jl")

# Define the frequency range
f = range(1e3, 10e3, length=100)
ω = 2π .* f

# Compute the sound velocities of the porous domain by using the attenuation coefficient
CP_attenuation = C_P.(ω)

# Compute the wavenumbers of the porous domain by using directly the wavenumber
wavenumber(ω) = compute_wave_properties_check(ω, sediment(ω))
aux = wavenumber.(ω)
kl_real = [r[1] for r in aux]
kl_imag = [r[2] for r in aux]

# Compute the sound velocities of the porous domain by using the wavenumbers
CP_wavenumber = ω./(kl_real - 1im * kl_imag)

# Compute the difference between the both sound velocities
diff = abs.(CP_wavenumber - CP_attenuation)./abs.(CP_wavenumber)

norm(diff)

@assert norm(diff) < 1e-2 