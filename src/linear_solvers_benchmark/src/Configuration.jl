"""
This file contains the configuration parameters for the simulation.
"""

# Load packages to use special functions in the calculation of the imaginary part of the bulk modulus
using SpecialFunctions

# Domain parametrization
L = 1. # Horizontal length of the domain [m]
t_P = 0.4 # Thickness of the porous domain [m]
t_F = 0.2 # Thickness of the fluid domain [m]
t = 0.25 # Length of the sonar transducer [m] # Is not being considered in the current implementation!!... with the translation of the solution this condition should not appear
h = 1 # Height of the water-column above t_F [m]
d_PML = 0.2 # Thickness of the PML layer [m]

# Sonar parametrization using the end-fire configuration
xᵦ = 0. # Horizontal coordinate of the sonar transducer [m] 
yᵦ = t_P + t_F + d_PML + h # Vertical coordinate of the sonar transducer [m]
Nₛ = 2 # Number of points to simulate the sonar transducer seen as Dirac's deltas [-]
Lₛ = 1*0.015 # Length of the sonar transducer [m]
Δy = Nₛ > 1 ? Lₛ/(Nₛ-1) : 0 # Separation between the points of the sonar discretization points [m]
Aₛ = 1/Nₛ # Amplitude of the sonar transducer [Pa]

# Scattering objects parametrization
    # Circles parametrization
    N_rigid_circles = 20 # Number of rigid circles in the porous domain [-]
    N_porous_circles = 0 # Number of porous circles in the porous domain [-]
    r_0 = 1.25e-2 # Interior radius of the circle [m]
    σ_r = 0.005 # Standard deviation of the interior radius of the circle [m]
    tol_circle = 0.5 * r_0 # Tolerance for the circles to avoid collisions with the boundaries of the physical domain
    
    # Clamps parametrization
    N_rigid_clamps = 0 # Number of rigid clamps in the porous domain [-]
    N_porous_clamps = 0 # Number of porous clamps in the porous domain [-]
    r_0i = 0.04 # Interior radius of the clamp [m]
    σ_ri = 0.05e-2 # Standard deviation of the interior radius of the clamp [m]
    r_0e = 0.042 # Exterior radius of the clamp [m] 
    σ_re = 0.05e-2 # Standard deviation of the interior radius of the clamp [m]
    θ_o_min = 0 # Minimum opening angle of the clamp [rad]
    θ_o_max = pi/9 # Maximum opening angle of the clamp [rad]
    θ_b_min = pi/9 # Minimum bisector angle of the clamp [rad]
    θ_b_max = pi # Maximum bisector angle of the clamp [rad]
    tol_clamp = 0.5 * r_0e # Tolerance for the clamps to avoid collisions with the boundaries of the physical domain
    
    # Ellipses parametrization
    N_rigid_ellipses = 20 # Number of rigid ellipses in the porous domain [-]
    a_0 = 0.03 # Major semi-axis of the ellipse [m]
    σ_a =  0.05e-2 # Standard deviation of the major semi-axis of the ellipse [m]
    b_0 = a_0/2 # Minor semi-axis of the ellipse [m]
    σ_b =  0.02e-2 # Standard deviation of the minor semi-axis of the ellipse [m]
    e_min = 0.002 # Minimum espesor of the ellipse [-]
    e_max = 0.002 # Maximum espesor of the ellipse [-]
    θ_el_min = pi/15 # Minimum angle of the ellipse [rad]
    θ_el_max = pi/15 # Maximum angle of the ellipse [rad]
    α_min = 0 # Minimum angle of rotation of the ellipse [rad]
    α_max = 2*pi # Maximum angle of rotation the ellipse [rad]
    tol_ellipse = 2 * a_0 # Tolerance for the ellipses to avoid collisions with the boundaries of the physical domain

# Frequency and angular frequency
f = 15e3 # Frequency [Hz]
ω = 2 * π * f # Angular frequency [rad/s]

# Transducer pressure [Pa]
P_0 = 1


# Physical properties of physical domains
    
    # Define the physical properties of the fluid domain, supposed to be uniform in the water column
    ρ_F(ω) = 1000. # Mass density of the fluid [kg/m^3]
    c_F(ω) = 1432. # Speed of sound in the fluid [m/s]
    
    # Define the physical properties of the porous domain
    ρᵣ = 2650. # Mass density of the grains of the porous domain [kg/m^3]
    α = 0.34 # Porosity of the porous domain [-]
    c = 1000. # Speed of sound in the porous domain [m/s]
    Kbᵣ = (1-α) * ρᵣ * c^2 # Real part of the bulk modulus of the porous domain [Pa]
    q = 0.074  # Log-decrement for longitudinal waves in the porous domain [-]
    Kbᵢ(ω) = Kbᵣ/(1 +
                  q*(sin(ω/c*(π/2 - sinint(ω/c)))-cos(ω/c*(cosint(ω/c)))) +
                  1im*q*(cos(ω/c*(sinint(ω/c)-π/2))-sin(ω/c*(cosint(ω/c))))
                  )  # Imaginary part of the bulk modulus of the porous domain [Pa]
    ρ_P(ω) = (1-α) * ρᵣ # Mass density of the porous domain [kg/m^3]
    # c_P(ω) = sqrt(K_P(ω)/ρ_P(ω)) # Speed of sound in the porous domain [m/s]
    c_P(ω) = 3500. # Speed of sound in the porous domain [m/s]

    # Define the physical properties of the circle scattering object when is considered porous
    ρ_S_circle(ω) = 3000. # Mass density of the circle scattering object [kg/m^3]
    c_S_circle(ω) = 500. # Speed of sound in the circle scattering object [m/s]  

    # Define the physical properties of the clamp scattering object when is considered porous
    ρ_S_clamp(ω) = 4000. # Mass density of the circle scattering object [kg/m^3]
    c_S_clamp(ω) = 3000. # Speed of sound in the circle scattering object [m/s]  

# Define the bulk modulus for the fluid, porous and scattering domains
K_F(ω) = ρ_F(ω)*c_F(ω)^2 # [Pa]
# K_P(ω) = Kbᵣ + Kbᵢ(ω) # [Pa]
K_P(ω) = ρ_P(ω)*c_P(ω) # [Pa]
K_S_clamp(ω) = ρ_S_clamp(ω) * c_S_clamp(ω)^2
K_S_circle(ω) = ρ_S_circle(ω) * c_S_circle(ω)^2

# Define the wavenumbers for the fluid, porous and scattering domains
k_F(ω) = ω/c_F(ω) # [rad/m]
k_P(ω) =  ω/c_P(ω) # [rad/m]
k_S_circle(ω) = ω/c_S_circle(ω)
k_S_clamp(ω) = ω/c_S_clamp(ω)




# PML parameters for the quadratic profile
R_PML = 1e-5 # Tolerence for PML reflection
σ_0 = -3/4*log(R_PML)/d_PML # PML coefficient


  


