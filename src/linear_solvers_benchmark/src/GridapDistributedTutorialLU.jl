using Gridap
using GridapGmsh
using GridapPETSc
using GridapDistributed
using PartitionedArrays

function main(ranks)

model = GmshDiscreteModel(ranks,"./data/demo.msh")
order = 1
dirichlet_tags = ["boundary1","boundary2"]
u_boundary1(x) = (0.0 + 1.0im)*x[2]
u_boundary2(x) = (1.0 + 0.0im)*x[1]
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,dirichlet_tags=dirichlet_tags, vector_type=Vector{ComplexF64})
U = TrialFESpace(V,[u_boundary1,u_boundary2])
Ω = Interior(model)
dΩ = Measure(Ω,2*order)
a(u,v) = (1.0+0.0im)*∫( ∇(u)⋅∇(v) )dΩ
l(v) = 0.0 + 0.0im
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
# writevtk(Ω, "results/demo", cellfields=["Real(uh)"=>real(uh), "Imag(uh)"=>imag(uh)])

end

with_mpi() do distribute 
  ranks = distribute_with_mpi(LinearIndices((4,)))
  main(ranks)
end