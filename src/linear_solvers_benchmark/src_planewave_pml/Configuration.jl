"""
This file contains the configuration parameters for the simulation.
"""

# Domain parametrization
L = 1 # Horizontal length of the domain [m]
t_P = 0.5 # Thickness of the porous domain [m]
d_PML = 0.2 # Thickness of the PML layer [m]

# Frequency and angular frequency
f_max = 10e3 # Frequency [Hz]
ω_max = 2 * π * f_max # Angular frequency [rad/s]

# Phyisical properties
ρ_P(ω) = 2000. # Mass density of the porous [kg/m^3]
c_P(ω) = 2000. # Speed of sound in the porous [m/s]

# Define the bulk modulus and wavenumber for the fluid and porous domains
K_P(ω) = ρ_P(ω) *c_P(ω)^2 # [Pa]
k_P(ω) =  ω/c_P(ω) # [rad/m]

# Transducer pressure [Pa]
P_0 = 5e5 + 5e5im 

# PML parameters for the quadratic profile
R_PML = 1e-5 # Tolerence for PML reflection
σ_0 = -3/4*log(R_PML)/d_PML # PML coefficient