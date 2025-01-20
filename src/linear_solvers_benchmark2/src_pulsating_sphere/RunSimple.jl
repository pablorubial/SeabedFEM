"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using  GridapGmsh

# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

labels = get_face_labeling(model)

Ω = Interior(model) # Computational domain
xp = get_physical_coordinate(Ω)
fun(x) = (x[1] > 0 && x[2]-H/2 > 0) ? 1000 : (x[1] < 0 && x[2]-H/2 < 0) ? 2000 : (x[1] > 0 && x[2]-H/2 < 0) ? 3000 : 4000

u = fun ∘ xp

writevtk(Ω,"./results/result_novel.vtu", cellfields=[ "u"=>u])

