using Gridap
using GridapGmsh
using GridapPETSc
using GridapDistributed
using PartitionedArrays



# Charge the model
model = GmshDiscreteModel("./data/demo.msh")

# Define the finite element space: R-T of order 1
order = 1
dirichlet_tags = ["boundary1"]
u_boundary1(x) = VectorValue(x[1], x[2], x[3])
reffe = ReferenceFE(raviart_thomas, Float64, order)
V = TestFESpace(model,reffe, conformity=:Hdiv, dirichlet_tags=dirichlet_tags, vector_type=Vector{ComplexF64})
U = TrialFESpace(V, [u_boundary1])

# Define the interior domain
Ω = Interior(model)
dΩ = Measure(Ω, 2*order)

# Define the bilinear form
a(u, v) = (1.0+0.0im)*∫((u⋅v)*6000^2)dΩ - ∫(divergence(u)*divergence(v))dΩ
l(v) = 0.0 + 0.0im

op = AffineFEOperator(a,l,U,V)
uh = solve(op)
writevtk(Ω, "results/demoRT", cellfields=["Real(uh)"=>real(uh), "Imag(uh)"=>imag(uh)])

