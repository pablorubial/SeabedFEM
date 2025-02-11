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
h1 = c/(15*f) # See how to extract this value from the value that Andres specified to me: h = (c/f)/15

# Corners of Porous Domain
p1 = gmsh.model.geo.addPoint(-L/2, 0, 0, h1)
p2 = gmsh.model.geo.addPoint(L/2, 0, 0, h1)
p3= gmsh.model.geo.addPoint(L/2, H, 0, h1) 
p4 = gmsh.model.geo.addPoint(-L/2, H, 0, h1)

# Corners of the Bottom PML
    # Bottom Left Part
    p5 = gmsh.model.geo.addPoint(-L/2-d_PML, 0, 0, h1)
    p6 = gmsh.model.geo.addPoint(-L/2-d_PML, -d_PML, 0, h1)
    # Bottom Central Part
    p7 = gmsh.model.geo.addPoint(-L/2, -d_PML, 0, h1)
    p8 = gmsh.model.geo.addPoint(L/2, -d_PML, 0, h1)
    # Bottom Right Part
    p9 = gmsh.model.geo.addPoint(L/2 + d_PML, 0, 0, h1)
    p10 = gmsh.model.geo.addPoint(L/2 + d_PML, -d_PML, 0, h1)


# Corners of the Top PML
    # Top Right Part
    p11 = gmsh.model.geo.addPoint(L/2+d_PML, H, 0, h1)
    p12 = gmsh.model.geo.addPoint(L/2 + d_PML, H+d_PML, 0, h1)
    # Top Central Part
    p13 = gmsh.model.geo.addPoint(L/2, H+d_PML, 0, h1)
    p14 = gmsh.model.geo.addPoint(-L/2, H+d_PML, 0, h1)
    # Top Left Part
    p15 = gmsh.model.geo.addPoint(-L/2-d_PML, H+d_PML, 0, h1)
    p16 = gmsh.model.geo.addPoint(-L/2-d_PML, H, 0, h1)

# Rock points
pr1 = gmsh.model.geo.addPoint(x_0 + r, y_0, 0, h1)
pr2 = gmsh.model.geo.addPoint(x_0, y_0, 0, h1)
pr3 = gmsh.model.geo.addPoint(x_0 - r, y_0, 0, h1)


## EDGES

# Edges of the Porous
l_P1 = gmsh.model.geo.add_line(p1, p2)
l_P2 = gmsh.model.geo.add_line(p2, p3)
l_P3 = gmsh.model.geo.add_line(p3, p4)
l_P4 = gmsh.model.geo.add_line(p4, p1)

# Edges of the Bottom PML
    # Left 
    l_pml_b1 = gmsh.model.geo.add_line(p1, p5)
    l_pml_b2 = gmsh.model.geo.add_line(p5, p6)
    l_pml_b3 = gmsh.model.geo.add_line(p6, p7)
    l_pml_b4 = gmsh.model.geo.add_line(p7, p1)
    # Center
    c_pml_b1 = -l_P1
    c_pml_b2 = -l_pml_b4
    c_pml_b3 = gmsh.model.geo.add_line(p7, p8)
    c_pml_b4 = gmsh.model.geo.add_line(p8, p2)    
    # Right Part
    r_pml_b1 = gmsh.model.geo.add_line(p9, p2)
    r_pml_b2 = -c_pml_b4
    r_pml_b3 = gmsh.model.geo.add_line(p8, p10)
    r_pml_b4 = gmsh.model.geo.add_line(p10, p9)

# Edges of the  Lateral PML
    # Left Part
    l_pml_l1 = gmsh.model.geo.add_line(p4, p16)
    l_pml_l2 = gmsh.model.geo.add_line(p16, p5)
    l_pml_l3 = -l_pml_b1
    l_pml_l4 = -l_P4 
    # Right Part
    r_pml_l1 = gmsh.model.geo.add_line(p11, p3)
    r_pml_l2 = -l_P2
    r_pml_l3 = -r_pml_b1
    r_pml_l4 = gmsh.model.geo.add_line(p9, p11)
# Edges of the top PML
    # Left Part
    l_pml_t1 = gmsh.model.geo.add_line(p14, p15)
    l_pml_t2 = gmsh.model.geo.add_line(p15, p16)
    l_pml_t3 = -l_pml_l1
    l_pml_t4 = gmsh.model.geo.add_line(p4, p14)
    # Center Part
    c_pml_t1 = gmsh.model.geo.add_line(p13, p14)
    c_pml_t2 = -l_pml_t4
    c_pml_t3 = -l_P3
    c_pml_t4 = gmsh.model.geo.add_line(p3, p13)
    # Right Part
    r_pml_t1 = gmsh.model.geo.add_line(p12, p13)
    r_pml_t2 = -c_pml_t4
    r_pml_t3 = -r_pml_l1
    r_pml_t4 = gmsh.model.geo.add_line(p11, p12)
# Edges of the rock
l_r1 = gmsh.model.geo.addCircleArc(pr1, pr2, pr3) # Rock_circle
l_r2 = gmsh.model.geo.addCircleArc(pr3, pr2, pr1) # Rock_circle

# Curve Loops
    # Physical Domain
    cl_P = gmsh.model.geo.add_curve_loop([l_P1, l_P2, l_P3, l_P4, -l_r1, -l_r2])
    # Bottom PML
        # Left Part
        cl_PML_lb = gmsh.model.geo.add_curve_loop([l_pml_b1, l_pml_b2, l_pml_b3, l_pml_b4])
        # Center Part
        cl_PML_cb = gmsh.model.geo.add_curve_loop([c_pml_b1, c_pml_b2, c_pml_b3, c_pml_b4])
        # Right Part
        cl_PML_rb = gmsh.model.geo.add_curve_loop([r_pml_b1, r_pml_b2, r_pml_b3, r_pml_b4])
    # Lateral PML
        # Left Part
        cl_PML_ll = gmsh.model.geo.add_curve_loop([l_pml_l1, l_pml_l2, l_pml_l3, l_pml_l4])
        # Right Part
        cl_PML_rl = gmsh.model.geo.add_curve_loop([r_pml_l1, r_pml_l2, r_pml_l3, r_pml_l4])
    # Top PML
        # Left Part
        cl_PML_lt = gmsh.model.geo.add_curve_loop([l_pml_t1, l_pml_t2, l_pml_t3, l_pml_t4])
        # Center Part
        cl_PML_ct = gmsh.model.geo.add_curve_loop([c_pml_t1, c_pml_t2, c_pml_t3, c_pml_t4])
        # Right Part
        cl_PML_rt = gmsh.model.geo.add_curve_loop([r_pml_t1, r_pml_t2, r_pml_t3, r_pml_t4])
# cl_rock = gmsh.model.geo.add_curve_loop([-l_r1, -l_r2])




# Surfaces
# s_P = gmsh.model.geo.addPlaneSurface([cl_P, cl_rock])
s_P = gmsh.model.geo.addPlaneSurface([cl_P])
s_PML_lb = gmsh.model.geo.addPlaneSurface([cl_PML_lb])
s_PML_cb = gmsh.model.geo.addPlaneSurface([cl_PML_cb])
s_PML_rb = gmsh.model.geo.addPlaneSurface([cl_PML_rb])
s_PML_ll = gmsh.model.geo.addPlaneSurface([cl_PML_ll])
s_PML_rl = gmsh.model.geo.addPlaneSurface([cl_PML_rl])
s_PML_lt = gmsh.model.geo.addPlaneSurface([cl_PML_lt])
s_PML_ct = gmsh.model.geo.addPlaneSurface([cl_PML_ct])
s_PML_rt = gmsh.model.geo.addPlaneSurface([cl_PML_rt])

s_list = [s_PML_lb, s_PML_cb, s_PML_rb, s_PML_ll, s_PML_rl, s_PML_lt, s_PML_ct, s_PML_rt]

# Synchronize model before meshing
gmsh.model.geo.synchronize()

#Set physical groups for 1D entities
l_rock = gmsh.model.addPhysicalGroup(1, [l_r1, l_r2])
l_pml = gmsh.model.addPhysicalGroup(1, [l_pml_b2, l_pml_b3, c_pml_b3, r_pml_b3, r_pml_b4, r_pml_l4, r_pml_t4, r_pml_t1, c_pml_t1, l_pml_t1, l_pml_t2, l_pml_l2])

# Set physical groups for 2D entities
f_P = gmsh.model.addPhysicalGroup(2, [s_P])
l_PML = gmsh.model.addPhysicalGroup(2, [s_PML_lb, s_PML_lt, s_PML_rb, s_PML_rt, r_pml_l4])
c_PML = gmsh.model.addPhysicalGroup(2, [s_PML_cb, s_PML_ct])
ll_pml = gmsh.model.addPhysicalGroup(2, [s_PML_ll, s_PML_rl])


# Set physical names for 1D entities
gmsh.model.setPhysicalName(1, l_rock, "source")
gmsh.model.setPhysicalName(1, l_pml, "sides")

# Set physical names for 2D entities
gmsh.model.setPhysicalName(2, f_P, "physicalDomain")
gmsh.model.setPhysicalName(2, l_PML, "cornerPML")
gmsh.model.setPhysicalName(2, c_PML, "centralPML")
gmsh.model.setPhysicalName(2, ll_pml, "lateralPML")


# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Write mesh to file
gmsh.write("./data/mesh.msh")
gmsh.finalize()

# # Convert the mesh to the json Gridap format
# model = DiscreteModelFromFile("./data/pml_mesh_quad.msh")
# # Write the mesh to a vtk file
# writevtk(model,"./results/pml_mesh")
# # Save the mesh in the json format
# fn = "./data/pml_mesh.json"
# to_json_file(model,fn)