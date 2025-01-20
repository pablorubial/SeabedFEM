"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap.Geometry
using GridapDistributed
using GridapDistributed.CellData
using PartitionedArrays
using Gridap
using Gridap.Fields
using GridapGmsh
using GridapPETSc

function main(parts)

    include("./Configuration.jl")

    model = GmshDiscreteModel(parts, "data/mesh.msh")

    Ω = Triangulation(model) # Computational domain
    # xp = get_physical_coordinate(Ω)
    # labels = get_face_labeling(model)
    # # fields = map_parts(labels.labels) do labels
    # #     dimension = 3
    # #     tags = get_face_tag(labels,dimension)
    # #     CellField(tags)
    # # end
    # # tags = DistributedCellField(tags)
    fun(x) = (x[1] > 0 && x[2]-H/2 > 0) ? 1000 : (x[1] < 0 && x[2]-H/2 < 0) ? 2000 : (x[1] > 0 && x[2]-H/2 < 0) ? 3000 : 4000
    
    # fields = map(local_views(Ω)) do trians
    #     GenericField(f)
    # end

    # labels = get_face_labeling(model)

    field = map(local_views(Ω)) do domains
        xp = get_physical_coordinate(domains)
    end

    # println(typeof(field))

    # u = fun ∘ field

    u = GridapDistributed.DistributedCellField(u)

    # u = CellField(fun, Ω)

    println()

    println(typeof(u))

    # fun = GridapDistributed.DistributedCellField(fields)

    # # u = fun ∘ xp

    # writevtk(Ω,"./results/result_simple", cellfields=[ "u"=>u])
    
end
    
with_mpi() do distribute 
    ranks = distribute_with_mpi(LinearIndices((4,)))
    main(ranks)
end