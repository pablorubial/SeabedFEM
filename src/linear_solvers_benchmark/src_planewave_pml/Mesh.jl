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
gmsh.model.add("domain_with_pml")

# Mesh size depending on the angular frequency
h = c_P(ω_max)/(20*f_max) # See how to extract this value from the value that Andres specified to me: h = (c/f)/15


# Corners of Porous Domain
p1 = gmsh.model.geo.addPoint(-L/2, 0, 0, h)
p2 = gmsh.model.geo.addPoint(L/2, 0, 0, h)
p3 = gmsh.model.geo.addPoint(L/2, t_P, 0, h) 
p4 = gmsh.model.geo.addPoint(-L/2, t_P, 0, h)

# Corners of the Bottom PML
p5 = gmsh.model.geo.addPoint(-L/2, -d_PML, 0, h)
p6 = gmsh.model.geo.addPoint(L/2, -d_PML, 0, h)


# Edges of the Porous
l_P1 = gmsh.model.geo.add_line(p1, p2)
l_P2 = gmsh.model.geo.add_line(p2, p3)
l_P3 = gmsh.model.geo.add_line(p3, p4)
l_P4 = gmsh.model.geo.add_line(p4, p1)


# Edges of the Bottom PML
l_PML1 = gmsh.model.geo.add_line(p1, p5)
l_PML2 = gmsh.model.geo.add_line(p5, p6)
l_PML3 = gmsh.model.geo.add_line(p6, p2)


# Curve Loops
cl_dom = gmsh.model.geo.add_curve_loop([l_P1, l_P2, l_P3, l_P4])
cl_pml = gmsh.model.geo.add_curve_loop([l_PML1, l_PML2, l_PML3, -l_P1])

# Surfaces
s_P = gmsh.model.geo.addPlaneSurface([cl_dom])
s_PML = gmsh.model.geo.addPlaneSurface([cl_pml])

# Synchronize model before meshing
gmsh.model.geo.synchronize()

# Set physical groups for 1D entities
f_transducer = gmsh.model.addPhysicalGroup(1, [l_P3])
f_dirichlet = gmsh.model.addPhysicalGroup(1, [l_P2, l_P4, l_PML1, l_PML2, l_PML3])

# Set physical groups for 2D entities
f_P = gmsh.model.addPhysicalGroup(2, [s_P])
f_PML = gmsh.model.addPhysicalGroup(2, [s_PML])


# Set physical names for 1D entities
gmsh.model.setPhysicalName(1, f_transducer, "transducer")
gmsh.model.setPhysicalName(1, f_dirichlet, "sides")


# Set physical names for 2D entities
gmsh.model.setPhysicalName(2, f_P, "physical_domain")
gmsh.model.setPhysicalName(2, f_PML, "PML")

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Write mesh to file
gmsh.write("./data/planewave.msh")
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/planewave.msh")
# Write the mesh to a vtk file
writevtk(model,"./results/planewave")