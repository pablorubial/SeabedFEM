"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh

# using GridapPETSc

# Define all the neccesary functions to compute the solution
# Define the tensors H, H^-1 and the Jacobian for the fluid PML only in the horizontal direction (quadratic profile)
include("./Configuration.jl")
include("./ExactSolution.jl")

u(x) = exact_solution(x)

# Load the mesh
model = GmshDiscreteModel("data/manufactured.msh")

# Define the finite element space: Lagrange of order 1
k = 1
reffe = ReferenceFE(lagrangian, Float64,k)
V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags=["left", "right", "top", "bottom"], vector_type=Vector{ComplexF64})
U = TrialFESpace(V, [u, u, u, u])

Ω = Triangulation(model) # Computational domain

# Define the measure for the fluid and porous domains
degree = 2

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

f(x) = -2*1im- 2

a(u,v) = ∫( (1.0+0.0im)*(∇(v)⊙∇(u)) )*dΩ
b(v) = ∫( v*f )*dΩ

op = AffineFEOperator(a,b,U,V)

uh = solve(op)

e = u - uh

# uh_real = CellField(V, real(uh.cell_dof_values))
# uh_imag = CellField(V, imag(uh.cell_dof_values))
u_ex = CellField(u, Ω)


writevtk(Ω,"./results/results_manufactured/results.vtu", cellfields=["Re(uh)"=>real(uh), "Im(u_h)"=>imag(uh),
                                                                     "Re(u)"=>real(u_ex), "Im(u)"=>imag(u_ex),
                                                                     "Re(error)"=>real(e), "Im(error)"=>imag(e)])  