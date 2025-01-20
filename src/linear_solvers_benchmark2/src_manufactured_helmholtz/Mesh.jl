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
gmsh.model.add("ManufacturedSolution")
# h = c/(20*freq) # See how to extract this value from the value that Andres specified to me: h = (c/f)/15
h = 0.005


# Corners of Porous Domain
p_P1 = gmsh.model.geo.addPoint(0, 0, 0, h)
p_P2 = gmsh.model.geo.addPoint(L, 0, 0, h)
p_P3 = gmsh.model.geo.addPoint(L, t_P, 0, h) 
p_P4 = gmsh.model.geo.addPoint(0, t_P, 0, h)

# Edges of the Porous
l_P1 = gmsh.model.geo.add_line(p_P1, p_P2)
l_P2 = gmsh.model.geo.add_line(p_P2, p_P3)
l_P3 = gmsh.model.geo.add_line(p_P3, p_P4)
l_P4 = gmsh.model.geo.add_line(p_P4, p_P1)
# Curve Loops
cl_P = gmsh.model.geo.add_curve_loop([l_P1, l_P2, l_P3, l_P4])

# Surfaces
s_P = gmsh.model.geo.addPlaneSurface([cl_P])

# Synchronize model before meshing
gmsh.model.geo.synchronize()

# Set physical groups for 1D entities
f_bottom = gmsh.model.addPhysicalGroup(1, [l_P1])
f_right = gmsh.model.addPhysicalGroup(1, [l_P2])
f_top = gmsh.model.addPhysicalGroup(1, [l_P3])
f_left = gmsh.model.addPhysicalGroup(1, [l_P4])

# Set physical groups for 2D entities
f_P = gmsh.model.addPhysicalGroup(2, [s_P])


# Set physical names for 1D entities
gmsh.model.setPhysicalName(1, f_bottom, "bottom")
gmsh.model.setPhysicalName(1, f_right, "right")
gmsh.model.setPhysicalName(1, f_top, "top")
gmsh.model.setPhysicalName(1, f_left, "left")

# Set physical names for 2D entities
gmsh.model.setPhysicalName(2, f_P, "domain")



# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Write mesh to file
gmsh.write("./data/manufactured_helmholtz.msh")
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/manufactured_helmholtz.msh")
# Write the mesh to a vtk file
writevtk(model,"./results/results_manufactured_helmholtz/mesh")