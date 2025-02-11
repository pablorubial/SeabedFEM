"""
This file generates the mesh of the PML domain using Gmsh 
The mesh includes the fluid domain, the porous domain and the PML and it is centered at the origin.
The mesh is saved in the folder data in the format .msh and .json
"""
# Load the packages
using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh

# Input of the physical parameters neccesary to calculate the mesh size
include("Configuration.jl")

# Initialize gmsh
gmsh.initialize()
gmsh.model.add("fluid_and_porous_PML")

# Mesh size depending on the angular frequency
h1 = c_F(ω)/(20*f) # See how to extract this value from the value that Andres specified to me: h = (c/f)/15
h2 = c_P(ω)/(20*f) # See how to extract this value from the value that Andres specified to me: h = (c/f)/15

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

# Corrners of the Transducer
p_T1 = gmsh.model.geo.addPoint(t/2, t_P+t_F, 0, h1)
p_T2 = gmsh.model.geo.addPoint(-t/2, t_P+t_F, 0, h1)

# Corners of the Fluid PML
    # Left Part
    #p_PML_F1 = gmsh.model.geo.addPoint(-L/2-d_PML, t_P, 0, h1)
    p_PML_F2 = gmsh.model.geo.addPoint(-L/2-d_PML, t_P+t_F, 0, h1)
    # Right Part
    #p_PML_F3 = gmsh.model.geo.addPoint(L/2+d_PML, t_P, 0, h1)
    p_PML_F4 = gmsh.model.geo.addPoint(L/2+d_PML, t_P+t_F, 0, h1)

# Corners of the Bottom PML
p_PML_B1 = gmsh.model.geo.addPoint(-L/2, -d_PML, 0, h2)
p_PML_B2 = gmsh.model.geo.addPoint(-L/2-d_PML, -d_PML, 0, h2)
#p_PML_B3 = gmsh.model.geo.addPoint(-L/2-d_PML, 0, 0, h2)
p_PML_B4 = gmsh.model.geo.addPoint(L/2, -d_PML, 0, h2)
p_PML_B5 = gmsh.model.geo.addPoint(L/2+d_PML, -d_PML, 0, h2)
#p_PML_B6 = gmsh.model.geo.addPoint(L/2+d_PML, 0, 0, h2)

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
l_T1 = gmsh.model.geo.add_line(p_F1, p_T1)
l_T2 = gmsh.model.geo.add_line(p_T1, p_T2)
l_F3 = gmsh.model.geo.add_line(p_T2, p_F2)
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


# Curve Loops
cl_P = gmsh.model.geo.add_curve_loop([l_P1, l_P2, l_P3, l_P4])
cl_PML_P1 = gmsh.model.geo.add_curve_loop([l_PML_P1, l_PML_P2, l_PML_P3, l_PML_P4])
cl_PML_P2 = gmsh.model.geo.add_curve_loop([l_PML_P5, l_PML_P6, l_PML_P7, l_PML_P8])
cl_F = gmsh.model.geo.add_curve_loop([l_F1, l_F2, l_T1, l_T2, l_F3, l_F4])
cl_PML_F1 = gmsh.model.geo.add_curve_loop([l_PML_F1, l_PML_F2, l_PML_F3, l_PML_F4])
cl_PML_F2 = gmsh.model.geo.add_curve_loop([l_PML_F5, l_PML_F6, l_PML_F7, l_PML_F8])
cl_PML_B1 = gmsh.model.geo.add_curve_loop([l_PML_B1, l_PML_B2, l_PML_B3, l_PML_B4])
cl_PML_B2 = gmsh.model.geo.add_curve_loop([l_PML_B5, l_PML_B6, l_PML_B7, l_PML_B8])
cl_PML_B3 = gmsh.model.geo.add_curve_loop([l_PML_B9, l_PML_B10, l_PML_B11, l_PML_B12])

# Surfaces
s_P = gmsh.model.geo.addPlaneSurface([cl_P])
s_PML_P1 = gmsh.model.geo.addPlaneSurface([cl_PML_P1])
s_PML_P2 = gmsh.model.geo.addPlaneSurface([cl_PML_P2])
s_F = gmsh.model.geo.addPlaneSurface([cl_F])
s_PML_F1 = gmsh.model.geo.addPlaneSurface([cl_PML_F1])
s_PML_F2 = gmsh.model.geo.addPlaneSurface([cl_PML_F2])
s_PML_B1 = gmsh.model.geo.addPlaneSurface([cl_PML_B1])
s_PML_B2 = gmsh.model.geo.addPlaneSurface([cl_PML_B2])
s_PML_B3 = gmsh.model.geo.addPlaneSurface([cl_PML_B3])

# Synchronize model before meshing
gmsh.model.geo.synchronize()

# Set physical groups for 1D entities
f_T = gmsh.model.addPhysicalGroup(1, [l_T1, l_T2, l_F3, l_PML_F8,l_PML_F2])
# f_T2 = gmsh.model.addPhysicalGroup(1, [l_PML_F8,l_PML_F2])
f_N = gmsh.model.addPhysicalGroup(1, [ l_PML_F3, l_PML_P3, l_PML_B4, l_PML_B10, l_PML_P7,  l_PML_F7])
f_B = gmsh.model.addPhysicalGroup(1, [l_PML_B1, l_PML_B5, l_PML_B9])

# Set physical groups for 2D entities
f_P = gmsh.model.addPhysicalGroup(2, [s_P])
f_PML_P = gmsh.model.addPhysicalGroup(2, [s_PML_P1, s_PML_P2])
f_F = gmsh.model.addPhysicalGroup(2, [s_F])
f_PML_F = gmsh.model.addPhysicalGroup(2, [s_PML_F1, s_PML_F2])
f_PML_B_s = gmsh.model.addPhysicalGroup(2, [s_PML_B1, s_PML_B3])
f_PML_B = gmsh.model.addPhysicalGroup(2, [s_PML_B2])
f_physical = gmsh.model.addPhysicalGroup(2, [s_P, s_F])

# Set physical names for 1D entities
gmsh.model.setPhysicalName(1, f_T, "transducer")
# gmsh.model.setPhysicalName(1, f_T2, "transducer_out")
gmsh.model.setPhysicalName(1, f_N, "sides")
gmsh.model.setPhysicalName(1, f_B, "bottom")

# Set physical names for 2D entities
gmsh.model.setPhysicalName(2, f_P, "porous_domain")
gmsh.model.setPhysicalName(2, f_PML_P, "porous_PML")
gmsh.model.setPhysicalName(2, f_F, "fluid_domain")
gmsh.model.setPhysicalName(2, f_PML_F, "fluid_PML")
gmsh.model.setPhysicalName(2, f_PML_B_s, "bottom_PML_lat")
gmsh.model.setPhysicalName(2, f_PML_B, "bottom_PML") 
gmsh.model.setPhysicalName(2, f_physical, "physical_domain")


# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Write mesh to file
gmsh.write("./data/planewave.msh")
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/planewave.msh")
# Write the mesh to a vtk file
writevtk(model,"./results/planewave_results/mesh")