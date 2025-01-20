"""
This file generates the mesh of the PML domain using Gmsh 
The mesh includes the fluid domain, the porous domain and the PML and it is centered at the origin.
The mesh is saved in the folder data in the format .msh and .json
"""
# Load the packages
using Random
Random.seed!(1234)
using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh
using Distributions
using Revise
# Input of the module that contain the matrix of circles object and the functions to generate it and check collisions between them
includet("./MarineFormsDict.jl")
using .MarineForm

# Input of the physical parameters neccesary to calculate the mesh size
include("Configuration.jl")

# Initialize gmsh
gmsh.initialize()
gmsh.model.add("fluid_and_porous_PML")

# Mesh size depending on the angular frequency
h1 = real(c_F(ω))/(15*f) 
h2 = real(c_P(ω))/(15*f) 
h_circle = real(c_S_circle(ω))/(20*f)
h_animal = h2/4

# Function to use polar coordinates
x(r, theta) =  r * cos(theta)
y(r, theta) = r * sin(theta)

# Corners of Porous Domain
p_P1 = gmsh.model.geo.addPoint(-L/2, 0, 0, h2)
p_P2 = gmsh.model.geo.addPoint(L/2, 0, 0, h2)
p_P3 = gmsh.model.geo.addPoint(L/2, t_P, 0, min(h1, h2)) 
p_P4 = gmsh.model.geo.addPoint(-L/2, t_P, 0, min(h1, h2))

# Corners of the Porous PML
    # Left Part
    p_PML_P1 = gmsh.model.geo.addPoint(-d_PML-L/2, 0, 0, h2)
    p_PML_P2 = gmsh.model.geo.addPoint(-d_PML-L/2, t_P, 0, min(h1, h2)) # == p_PML_F1
    # Right Part
    p_PML_P3 = gmsh.model.geo.addPoint(L/2+d_PML, 0, 0, h2)
    p_PML_P4 = gmsh.model.geo.addPoint(L/2+d_PML, t_P, 0, min(h1, h2)) # == p_PML_F3 

# Corrners of Fluid Domain
p_F1 = gmsh.model.geo.addPoint(L/2, t_P+t_F, 0, h1) 
p_F2 = gmsh.model.geo.addPoint(-L/2, t_P+t_F, 0, h1)

# Corners of the Fluid PML
    # Left Part
    p_PML_F2 = gmsh.model.geo.addPoint(-L/2-d_PML, t_P+t_F, 0, h1)
    # Right Part
    p_PML_F4 = gmsh.model.geo.addPoint(L/2+d_PML, t_P+t_F, 0, h1)

# Corners of the Bottom PML
p_PML_B1 = gmsh.model.geo.addPoint(-L/2, -d_PML, 0, h2)
p_PML_B2 = gmsh.model.geo.addPoint(-L/2-d_PML, -d_PML, 0, h2)
p_PML_B4 = gmsh.model.geo.addPoint(L/2, -d_PML, 0, h2)
p_PML_B5 = gmsh.model.geo.addPoint(L/2+d_PML, -d_PML, 0, h2)

# Corners of the Top PML
p_PML_T1 = gmsh.model.geo.addPoint(-L/2-d_PML, t_P+t_F+d_PML, 0, h1)
p_PML_T2 = gmsh.model.geo.addPoint(-L/2, t_P+t_F+d_PML, 0, h1)
p_PML_T3 = gmsh.model.geo.addPoint(L/2, t_P+t_F+d_PML, 0, h1)
p_PML_T4 = gmsh.model.geo.addPoint(L/2+d_PML, t_P+t_F+d_PML, 0, h1)


# Edges of the Porous
l_P1 = gmsh.model.geo.add_line(p_P1, p_P2)
l_P2 = gmsh.model.geo.add_line(p_P2, p_P3)
l_P3 = gmsh.model.geo.add_line(p_P3, p_P4)
l_P4 = gmsh.model.geo.add_line(p_P4, p_P1)

# Edges of the Porous PML
    # Left Part
    l_PML_P1 = -l_P4
    l_PML_P2 = gmsh.model.geo.add_line(p_P4, p_PML_P2)
    l_PML_P3 = gmsh.model.geo.add_line(p_PML_P2, p_PML_P1)
    l_PML_P4 = gmsh.model.geo.add_line(p_PML_P1, p_P1)
    # Right Part
    l_PML_P5 = -l_P2
    l_PML_P6 = gmsh.model.geo.add_line(p_P2, p_PML_P3)
    l_PML_P7 = gmsh.model.geo.add_line(p_PML_P3, p_PML_P4)
    l_PML_P8 = gmsh.model.geo.add_line(p_PML_P4, p_P3)

# Edges of the Fluid
l_F1 = -l_P3
l_F2 = gmsh.model.geo.add_line(p_P3, p_F1)
l_F3 = gmsh.model.geo.add_line(p_F1, p_F2)
l_F4 = gmsh.model.geo.add_line(p_F2, p_P4)

# Edges of the Fluid PML
    # Left Part
    l_PML_F1 = -l_F4
    l_PML_F2 = gmsh.model.geo.add_line(p_F2, p_PML_F2)
    l_PML_F3 = gmsh.model.geo.add_line(p_PML_F2, p_PML_P2)
    l_PML_F4 = -l_PML_P2

    # Right Part
    l_PML_F5 = -l_F2
    l_PML_F6 = -l_PML_P8
    l_PML_F7 = gmsh.model.geo.add_line(p_PML_P4, p_PML_F4)
    l_PML_F8 = gmsh.model.geo.add_line(p_PML_F4, p_F1)

# Edges of the Bottom PML
    # Left Part
    l_PML_B1 = gmsh.model.geo.add_line(p_PML_B2, p_PML_B1)
    l_PML_B2 = gmsh.model.geo.add_line(p_PML_B1, p_P1)
    l_PML_B3 = -l_PML_P4
    l_PML_B4 = gmsh.model.geo.add_line(p_PML_P1, p_PML_B2)
    # Center Part
    l_PML_B5 = gmsh.model.geo.add_line(p_PML_B1, p_PML_B4)
    l_PML_B6 = gmsh.model.geo.add_line(p_PML_B4, p_P2)
    l_PML_B7 = -l_P1
    l_PML_B8 = -l_PML_B2 # gmsh.model.geo.add_line(p_P1, p_PML_B1)
    # Right Part
    l_PML_B9 = gmsh.model.geo.add_line(p_PML_B4, p_PML_B5)
    l_PML_B10 = gmsh.model.geo.add_line(p_PML_B5, p_PML_P3)
    l_PML_B11 = -l_PML_P6
    l_PML_B12 = -l_PML_B6 # gmsh.model.geo.add_line(p_P2, p_PML_B4)

# Edges of the Top PML
    # Left Part
    l_PML_T1 = gmsh.model.geo.add_line(p_F2, p_PML_T2)
    l_PML_T2 = gmsh.model.geo.add_line(p_PML_T2, p_PML_T1)
    l_PML_T3 = gmsh.model.geo.add_line(p_PML_T1, p_PML_F2)
    l_PML_T4 = -l_PML_F2
    # Center Part
    l_PML_T5 = gmsh.model.geo.add_line(p_F1, p_PML_T3)
    l_PML_T6 = gmsh.model.geo.add_line(p_PML_T3, p_PML_T2)
    l_PML_T7 = -l_PML_T1
    l_PML_T8 = -l_F3 # gmsh.model.geo.add_line(p_P1, p_PML_B1)
    # Right Part
    l_PML_T9 = gmsh.model.geo.add_line(p_PML_F4, p_PML_T4)
    l_PML_T10 = gmsh.model.geo.add_line(p_PML_T4, p_PML_T3)
    l_PML_T11 = -l_PML_T5
    l_PML_T12 = -l_PML_F8 # gmsh.model.geo.add_line(p_P2, p_PML_B4)

# Curve Loops
cl_P = gmsh.model.geo.add_curve_loop([l_P1, l_P2, l_P3, l_P4])
cl_PML_P1 = gmsh.model.geo.add_curve_loop([l_PML_P1, l_PML_P2, l_PML_P3, l_PML_P4])
cl_PML_P2 = gmsh.model.geo.add_curve_loop([l_PML_P5, l_PML_P6, l_PML_P7, l_PML_P8])
cl_F = gmsh.model.geo.add_curve_loop([l_F1, l_F2, l_F3, l_F4])
cl_PML_F1 = gmsh.model.geo.add_curve_loop([l_PML_F1, l_PML_F2, l_PML_F3, l_PML_F4])
cl_PML_F2 = gmsh.model.geo.add_curve_loop([l_PML_F5, l_PML_F6, l_PML_F7, l_PML_F8])
cl_PML_B1 = gmsh.model.geo.add_curve_loop([l_PML_B1, l_PML_B2, l_PML_B3, l_PML_B4])
cl_PML_B2 = gmsh.model.geo.add_curve_loop([l_PML_B5, l_PML_B6, l_PML_B7, l_PML_B8])
cl_PML_B3 = gmsh.model.geo.add_curve_loop([l_PML_B9, l_PML_B10, l_PML_B11, l_PML_B12])
cl_PML_T1 = gmsh.model.geo.add_curve_loop([l_PML_T1, l_PML_T2, l_PML_T3, l_PML_T4])
cl_PML_T2 = gmsh.model.geo.add_curve_loop([l_PML_T5, l_PML_T6, l_PML_T7, l_PML_T8])
cl_PML_T3 = gmsh.model.geo.add_curve_loop([l_PML_T9, l_PML_T10, l_PML_T11, l_PML_T12])

# Definition of the objects of the bottom
definition = [
    (RigidCircle, Dict(
    :N => N_rigid_circles,
    :r_distribution => Normal(r_0, σ_r),
    :x_range => (-L/2 + (r_0 + 4*σ_r + tol_circle), L/2 - (r_0 + 4*σ_r + tol_circle)),
    :y_range => (0+(r_0 + 4*σ_r + + tol_circle), t_P-(r_0 + 4*σ_r + + tol_circle)))),
    
    (PorousCircle, Dict(
    :N => N_porous_circles,
    :r_distribution => Normal(r_0, σ_r),
    :x_range => (-L/2 + (r_0 + 4*σ_r + tol_circle), L/2 - (r_0 + 4*σ_r + tol_circle)),
    :y_range => (0+(r_0 + 4*σ_r + + tol_circle), t_P-(r_0 + 4*σ_r + + tol_circle)))),

    (RigidClamp, Dict(
    :N => N_rigid_clamps,
    :ri_distribution => Normal(r_0i, σ_ri),
    :re_distribution => Normal(r_0e, σ_re),
    :x_range => (-L/2 + (r_0e + 4*σ_re + tol_clamp), L/2 - (r_0e + 4*σ_re + tol_clamp)),
    :y_range => (0+(r_0e + 4*σ_re + tol_clamp), t_P-(r_0e + 4*σ_re + tol_clamp)),
    :θo_range => (θ_o_min, θ_o_max),
    :θb_range => (θ_b_min, θ_b_max))),
    
    (PorousClamp, Dict(
    :N => N_porous_clamps,
    :ri_distribution => Normal(r_0i, σ_ri),
    :re_distribution => Normal(r_0e, σ_re),
    :x_range => (-L/2 + (r_0e + 4*σ_re + tol_clamp), L/2 - (r_0e + 4*σ_re + tol_clamp)),
    :y_range => (0+(r_0e + 4*σ_re + tol_clamp), t_P-(r_0e + 4*σ_re + tol_clamp)),
    :θo_range => (θ_o_min, θ_o_max),
    :θb_range => (θ_b_min, θ_b_max))),

    (RigidEllipse, Dict(
        :N => N_rigid_ellipses,
        :a_distribution => Normal(a_0, σ_a),
        :b_distribution => Normal(b_0, σ_b),
        :e_range => (e_min, e_max),
        :x_range => (-L/2 + (a_0 + 4*σ_a + tol_ellipse), L/2 - (a_0 + 4*σ_a + tol_ellipse)),
        :y_range => (0+(a_0 + 4*σ_a + tol_ellipse), t_P-(a_0 + 4*σ_a + tol_ellipse)),
        :θo_range => (θ_el_min, θ_el_max),
        :α_range => (α_min, α_max))),
]

# Generate the animals
animals = generate_animals(definition)

rigid_circles = [animal for animal in animals if isa(animal, RigidCircle)]
porous_circles = [animal for animal in animals if isa(animal, PorousCircle)]
rigid_clamps = [animal for animal in animals if isa(animal, RigidClamp)]
porous_clamps = [animal for animal in animals if isa(animal, PorousClamp)]
rigid_ellipses = [animal for animal in animals if isa(animal, RigidEllipse)]

s_rigid_clamps = [generate_geometry(i, h_animal) for i in rigid_clamps]
s_porous_clamps = [generate_geometry(i, h_animal) for i in porous_clamps]
s_rigid_circles = [generate_geometry(i, h_animal) for i in rigid_circles]
s_porous_circles = [generate_geometry(i, h_circle) for i in porous_circles]
s_rigid_ellipses = [generate_geometry(i, h_animal) for i in rigid_ellipses]

cl_rigid_clamps = [i[1] for i in s_rigid_clamps]
cl_porous_clamps = [i[1] for i in s_porous_clamps]
cl_rigid_circles = [i[1] for i in s_rigid_circles]
cl_porous_circles = [i[1] for i in s_porous_circles]
cl_rigid_ellipses = [i[1] for i in s_rigid_ellipses]

l_rigid_clamps = [i[2] for i in s_rigid_clamps]
l_rigid_circles = [i[2] for i in s_rigid_circles]
l_rigid_ellipses = [i[2] for i in s_rigid_ellipses]
l_porous_circles = [i[2] for i in s_porous_circles]


# Surfaces
s_P = gmsh.model.geo.addPlaneSurface(vcat(cl_P, cl_rigid_clamps, cl_porous_clamps, cl_rigid_circles, cl_porous_circles, cl_rigid_ellipses))
s_circles = [gmsh.model.geo.addPlaneSurface([i]) for i in cl_porous_circles]
s_clamps = [gmsh.model.geo.addPlaneSurface([i]) for i in cl_porous_clamps]
s_PML_P1 = gmsh.model.geo.addPlaneSurface([cl_PML_P1])
s_PML_P2 = gmsh.model.geo.addPlaneSurface([cl_PML_P2])
s_F = gmsh.model.geo.addPlaneSurface([cl_F])
s_PML_F1 = gmsh.model.geo.addPlaneSurface([cl_PML_F1])
s_PML_F2 = gmsh.model.geo.addPlaneSurface([cl_PML_F2])
s_PML_B1 = gmsh.model.geo.addPlaneSurface([cl_PML_B1])
s_PML_B2 = gmsh.model.geo.addPlaneSurface([cl_PML_B2])
s_PML_B3 = gmsh.model.geo.addPlaneSurface([cl_PML_B3])
s_PML_T1 = gmsh.model.geo.addPlaneSurface([cl_PML_T1])
s_PML_T2 = gmsh.model.geo.addPlaneSurface([cl_PML_T2])
s_PML_T3 = gmsh.model.geo.addPlaneSurface([cl_PML_T3])

# Synchronize model before meshing
gmsh.model.geo.synchronize()

# Set physical groups for 1D entities
f_T = gmsh.model.addPhysicalGroup(1, [l_PML_F8, l_F3, l_PML_F2])
f_interface = gmsh.model.addPhysicalGroup(1, [l_F1, l_PML_P2, l_PML_P8])
f_N = gmsh.model.addPhysicalGroup(1, [ l_PML_F3, l_PML_P3, l_PML_B4, l_PML_B10, l_PML_P7,  l_PML_F7, l_PML_T9, l_PML_T10, l_PML_T6, l_PML_T2, l_PML_T3])
f_B = gmsh.model.addPhysicalGroup(1, [l_PML_B1, l_PML_B5, l_PML_B9])
f_animals = gmsh.model.addPhysicalGroup(1, vcat(l_rigid_circles..., l_rigid_clamps..., l_rigid_ellipses...))
f_porous_circles = gmsh.model.addPhysicalGroup(1, vcat(l_porous_circles...))


# Set physical groups for 2D entities
f_P = gmsh.model.addPhysicalGroup(2, [s_P, s_PML_B2, s_PML_B1, s_PML_B3, s_PML_P1, s_PML_P2])
f_F = gmsh.model.addPhysicalGroup(2, [s_F, s_PML_T1, s_PML_T3, s_PML_T2, s_PML_F1, s_PML_F2])
f_circles = gmsh.model.addPhysicalGroup(2, s_circles)
f_clamps = gmsh.model.addPhysicalGroup(2, s_clamps)
f_P_pysical = gmsh.model.addPhysicalGroup(2, [s_P])
f_F_physical = gmsh.model.addPhysicalGroup(2, [s_F])

# Set physical names for 1D entities
gmsh.model.setPhysicalName(1, f_N, "sides")
gmsh.model.setPhysicalName(1, f_B, "bottom")
gmsh.model.setPhysicalName(1, f_animals, "objects")
gmsh.model.setPhysicalName(1, f_interface, "interface")
gmsh.model.setPhysicalName(1, f_porous_circles, "circles_boundary")

# Set physical names for 2D entities
gmsh.model.setPhysicalName(2, f_P, "porous_domain")
gmsh.model.setPhysicalName(2, f_F, "fluid_domain")
gmsh.model.setPhysicalName(2, f_P_pysical, "porous_physical_domain")
gmsh.model.setPhysicalName(2, f_F_physical, "fluid_physical_domain")
gmsh.model.setPhysicalName(2, f_circles, "circles")
gmsh.model.setPhysicalName(2, f_clamps, "clamps") 

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Write mesh to file
gmsh.write("./data/mesh.msh")
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/mesh.msh")
# Write the mesh to a vtk file
writevtk(model,"./results/mesh")
# Save the mesh in the json format
# fn = "./data/pml_mesh.json"
# to_json_file(model,fn)
