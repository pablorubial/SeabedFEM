# Warning! Once the code is functional, is recommendable to define all the parameters with a 
# const value (i.e const a = 3), because this will optimize the use of Julia JIT

# Domain parametrization
L = 2 # Length of the domain [m]
H = 1 # Height of the domain [m]
d_PML = 0.15 # Thickness of the PML [m]
r = 0.1 # Radius of the rock [m]
x_0 = 0 # x-coordinate of the center of the source [m]
y_0 = H / 2 # y-coordinate of the center of the source [m]

# Check possible Domain
# if t > L
#     error("The length of the sonar cannot be bigger than the length of the domain")
# end

# Frequency parameters
f = 15e3
# ω = 2 * π * f

# Domains properties
ρ_F(ω) = 2000.
c_F(ω) = 8000.

# Transducer pressure
P_0 = 5e5im

# Wavenumber
k_F(ω) = ω / c_F(ω)

RPML = 1e-6 # Reflection coefficient of the PML
σ_0 = -3/4 * log(RPML)/d_PML
