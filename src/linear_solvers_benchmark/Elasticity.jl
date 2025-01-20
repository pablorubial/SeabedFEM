using Gridap.Geometry
using GridapDistributed
using GridapDistributed.CellData
using PartitionedArrays
using Gridap
using GridapGmsh

function main(parts)
      model = GmshDiscreteModel(parts,"/path/to/mesh")
      order = 1

      reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
      V0 = TestFESpace(model,reffe; conformity=:H1,dirichlet_tags=["surface_1","surface_2"], dirichlet_masks=[(true,false,false), (true,true,true)])

      g1(x) = VectorValue(0.005,0.0,0.0)
      g2(x) = VectorValue(0.0,0.0,0.0)

      # From functions `g1` and `g2`, we define the trial space as follows:

      U = TrialFESpace(V0,[g1,g2])

      degree = 2*order
      Ω = Triangulation(model)
      dΩ = Measure(Ω,degree)

      labels = get_face_labeling(model)
      fields = map_parts(labels.labels) do labels
             dimension = 3
             tags = get_face_tag(labels,dimension)
             CellField(tags,Ω)
         end
      tags = GridapDistributed.DistributedCellField(fields)
      alu_tag = get_tag_from_name(labels.labels.part,"material_1")

      function lame_parameters(E,ν)
           λ = (E*ν)/((1+ν)*(1-2*ν))
           μ = E/(2*(1+ν))
          (λ, μ)
      end

      E_alu = 70.0e9
      ν_alu = 0.33
      (λ_alu,μ_alu) = lame_parameters(E_alu,ν_alu)

      E_steel = 200.0e9
      ν_steel = 0.33
      (λ_steel,μ_steel) = lame_parameters(E_steel,ν_steel)

# Then, we define the function containing the constitutive law:

      function σ_bimat(ε,tag)
         if tag == alu_tag
             return λ_alu*tr(ε)*one(ε) + 2*μ_alu*ε
         else
              return λ_steel*tr(ε)*one(ε) + 2*μ_steel*ε
         end
      end

      a(u,v) = ∫( ε(v) ⊙ (σ_bimat∘(ε(u),tags) ))*dΩ 
      l(v) = 0

       op = AffineFEOperator(a,l,U,V0)
       uh = solve(op)

      #  writevtk(Ω,"demo_bimat",cellfields=["uh"=>uh,"epsi"=>ε(uh),"sigma"=>σ_bimat∘(ε(uh),tags)])
      writevtk(Ω,"demo_bimat",cellfields=["uh"=>uh,"epsi"=>ε(uh)])
end
partition = (1,3,2)
prun(main, mpi, partition)