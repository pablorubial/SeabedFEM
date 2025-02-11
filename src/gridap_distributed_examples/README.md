# GridapDistributed + GridapPETSc tests

This repository contains several examples of how to setup a `Gridap` case to be run in parallel by using its parallelization capabilities through `GridapDistributed` that performs all the process of mesh reading, assembly in parallel and `GridapPETSc` which solves the resultant linear system in parallel, by the used of `MPI`. It is recommended to look to the README file of the `../linear_solvers_benchmark` repository and follow the instructions to install and compile PETSc locally and link it with `GridapPETSc`.

There are three examples in three different folders:
* `src_manufactured_helmholtz`: which solves the Helmholtz equation by using the manufacturing solution methods of a given analytical solution. The corresponding running file is `RunDistribitedPETSc.jl`
* `src_pulsating_sphere`: which solves the problem `../benchmark2` in parallel. The corresponding running file is `RunDistribitedPETSc.jl`
* `src_planewave`: which solves the problem `../benchmark3` in parallel. The corresponding running file is `RunDistribitedPETSc.jl`