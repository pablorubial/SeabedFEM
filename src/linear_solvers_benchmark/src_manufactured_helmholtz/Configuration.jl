"""
This file contains the configuration parameters for the simulation.
"""

# Domain parametrization
L = 1 # Horizontal length of the domain [m]
t_P = 0.5 # Thickness of the porous domain [m]
# k_x = 5.75
# k_y = 3.45
c = 1000



# The maximum modes chosen are n_fr = 4 and m_fr = 4, this properties are used to mesh the domain
n_fr = 4
m_fr = 4

ω_r = c*π*sqrt((n_fr/L)^2 + (m_fr/t_P)^2)
k_r = ω_r/c
f_r = ω_r/(2π)

