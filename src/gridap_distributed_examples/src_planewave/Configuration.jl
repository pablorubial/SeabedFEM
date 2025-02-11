"""
This file contains the configuration parameters for the simulation.
"""

# Domain parametrization
L = 1 # Horizontal length of the domain [m]
t_P = 0.5 # Thickness of the porous domain [m]
t_F = 1.5 # Thickness of the fluid domain [m]
t = 0.25 # Length of the sonar transducer [m]
d_PML = 0.45 # Thickness of the PML layer [m]

# Frequency and angular frequency
f = 15e3 # Frequency [Hz]
ω = 2 * π * f # Angular frequency [rad/s]

# Phyisical properties
ρ_F(ω) = 1000. # Mass density of the fluid [kg/m^3]
c_F(ω) = 8000. # Speed of sound in the fluid [m/s]
ρ_P(ω) = 2000. # Mass density of the porous [kg/m^3]
c_P(ω) = 6000. # Speed of sound in the porous [m/s]

# Define the bulk modulus and wavenumber for the fluid and porous domains
K_F(ω) = ρ_F(ω)*c_F(ω)^2 # [Pa]
K_P(ω) = ρ_P(ω) *c_P(ω)^2 # [Pa]
k_F(ω) = ω/c_F(ω) # [rad/m]
k_P(ω) =  ω/c_P(ω) # [rad/m]

# Transducer pressure [Pa]
P_0 = 5e5 + 5e5im 

# PML parameters for the quadratic profile
R_PML = 1e-5 # Tolerence for PML reflection
σ_0 = -3/4*log(R_PML)/d_PML # PML coefficient