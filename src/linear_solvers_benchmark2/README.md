# Linking a PETSc local installation with GridapPETSc

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)]()



## Installation

This document will serve as a guide to try to link `GridapPETSc` installation with a local `PETSc` installation. Following the `GridapPETSc` documentation, it requires the `PETSC` library ([Portable, Extensible Toolkit for Scientific Computation](https://www.mcs.anl.gov/petsc/)) and `MPI` to work correctly. 

 - The first step is to clone the oficial PETSc repository in a local folder use to store the local installed software.
    ```bash
    mkdir $HOME/software
    cd $HOME/software
    git clone -b release https://gitlab.com/petsc/petsc.git petsc
    git pull # obtain new release fixes (since a prior clone or pull)
    ```
    Once the PETSc library have been cloned a folder called `petsc` should be visible inside the `$HOME/software` folder. The next step is compile it, to do this go into the `petsc` folder by doing:
    ```bash
    cd petsc
    ```

- There are not so much information about how to proceed to compile the library, some of the flags used has been obtained from the following documentation of an old Julia library [MOOSE.jl]([MOOSE.jl](https://friedmud.github.io/MOOSE.jl/installation/)). Many different compilations of PETSc can be live together in the `/$HOME /software/petsc` folder. For example, one can compile PETSc with a specific quantity of libraries or even more, customized to work with complex or real numbers. The `--PETSC_ARCH` flag establish the name of the folder where the specific compilated version of PETSc will be allocated. For example if one want to compile PETSc to work with real arithmetic, one flags that have worked in my case, are the following ones:

    ```bash
    # Configure PETSc
    python3 ./configure \ 
        --PETSC_ARCH=linux-gnu-real-64 \    
        --with-cc=mpicc \    
        --with-cxx=mpicxx \    
        --with-fc=mpif90 \  
        --with-scalar-type=real \  
        --with-debugging=no \    
        --with-pic=1 \    
        --with-shared-libraries=1 \    
        --download-metis=1 \    
        --download-parmetis=1 \    
        --download-superlu_dist=1 \    
        --download-mumps=1 \    
        --download-scalapack=1 \    
        --download-hypre=1 \    
        --download-cmake
    ```
    In the next step, the compilation usen make should be performed
    ```bash
    # Compile PETSc
    make PETSC_DIR=/home/pablo/software/petsc PETSC_ARCH=linux-gnu-real-64 all
    ```
    and finally the checking of the installed libraries to see if the installlation was succesfful
    ```bash
    # Check PETSc build
    make PETSC_DIR=/home/pablo/software/petsc PETSC_ARCH=linux-gnu-real-64 check
    ```
    At the same time, if one desires work with complex arithmetic it should be specified by the use of the flag `--with-scalar-type=complex` as follows:
    ```bash
    # Configure PETSc
    python3 ./configure \ 
        --PETSC_ARCH=linux-gnu-complex-64 \
        --with-cc=mpicc \
        --with-cxx=mpicxx \
        --with-fc=mpif90 \
        --with-scalar-type=complex \
        --with-pic=1 \
        --with-shared-libraries=1 \
        --with-debugging=no \
        --download-fblaslapack=1 \
        --download-metis=1 \
        --download-parmetis=1 \
        --download-superlu_dist=1 \
        --download-mumps=1 \
        --download-scalapack=1 \
        --download-hypre=1 \
        --download-elemental=1 \
        --download-cmake
    ```
    In the next step, the compilation usen make should be performed
    ```bash        
    # Compile PETSc
    make PETSC_DIR=/home/pablo/software/petsc PETSC_ARCH=linux-gnu-complex-64 all
    ```
    and finally the checking of the isntalled libraries to see if the installlation was succesfful
    ```bash
    # Check PETSc build
    make PETSC_DIR=/home/pablo/software/petsc PETSC_ARCH=linux-gnu-complex-64 check
    ```

    Once this two versions of PETSc have been compiled, inside the `petsc` library should be visible two folders with the name specified in the `--PETSC_ARCH` flag, in this case two folders with the names:

    ```
    linux-gnu-real-64 
    linux-gnu-complex-64 
    ```
    and if differents versions are compiled, they should ve visible inside the `$HOME  /software/petsc` folder. Following this steps the PETSc library is downloaded and compiled.
    
* **The chosen `PETSc` library needs to be configured with the `MPI` installation that the `MPI.jl` library is going to use**. In our case, we are telling to PETSc during the compilation to use the local MPI installation in the system with the flags `-with-cc, with-cxx, with-fc`. In Ubuntu 22.04 the default MPI installation of the system is OpenMPI. So when using the MPI Julia library, is neccesary to link the library with the system MPI library instead with the pre-built binaries that the by defaulf `MPI.jl` installation uses. These steps are quite good explained in the [MPI](https://juliaparallel.org/MPI.jl/stable/configuration/) official documentation. Following the instructions:
    ```bash
    # Create a Julia environment where the project is going to be developed
    julia --project=.
    ```
    Now, install the `MPIPreferences.jl` package that let to customize the MPI version used by the `MPI.jl` version.
    ```julia
    using Pkg
    Pkg.add("MPIPreferences")
    ```
    And now execute the following commands:
    ```julia
    using MPIPreferences
    MPIPreferences.use_system_binary()
    ```
    that will generate in the project folder a file with the name `LocalPreferences.toml` with the following content:
    ```bash
    [MPIPreferences]
    __clear__ = ["preloads_env_switch"]
    _format = "1.0"
    abi = "OpenMPI"
    binary = "system"
    cclibs = []
    libmpi = "/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so"
    mpiexec = "/usr/bin/mpiexec"
    preloads = []
    ```
    the only thing that the user should provide here is the path inside the system where the `libmpi` and `mpiexec` library and commands are located in the system. That in Ubuntu 22.04 seems to be that ones, since the probes done have been successful. With this, the user is specifying to the `MPI.jl` where the local MPI libray is located and also forcing to use it. So the `MPI.jl` library can be installed:
    ```julia
    Pkg.add("MPI")
    ```

- The next step is to install the `GridapPETSc.jl` library. By default, when the library is added to the project, it would use the PETSc Julia precompiled binaries that are available inside the `PETSc.jl` library. To enforce to `GridapPETSc.jl` to use our local PETSc installation we have to create an environment variable called `JULIA_PETSC_LIBRARY` and specify on it the PATH where the specific PETSc compiled library is located. For example, if we want to work with real arithmetic we have to add the following variable environments to our `~/.bashrc` file:
    ```bash
    export PETSC_DIR=/$HOME/software/petsc
    export PETSC_ARCH=linux-gnu-real-64
    export JULIA_PETSC_LIBRARY=/home/pablo/software/petsc/linux-gnu-real-64/lib/libpetsc.so
    ```
    Or if we want to work with complex arithmetics, changing the environment variables in the following way:
    ```bash
    export PETSC_DIR=/$HOME/software/petsc
    export PETSC_ARCH=linux-gnu-complex-64
    export JULIA_PETSC_LIBRARY=/home/pablo/software/petsc/linux-gnu-complex-64/lib/libpetsc.so
    ```
    With all the settings ready, the next step is to install the `GridapPETSc.jl` libray and build it to force it to use the local `PETSc` installation:
    ```julia
    Pkg.add("GridapPETSc")
    Pkg.build()
    ```

    With all these steps the user is able to use `GridapPETSc.jl` linked together with his own local `PETSc` installation:
    ```julia
    using GridapPETSc
    ```

## Notes of developers 

* All this proccedure have been tested for a native Ubuntu 22.04 installation, there is not any guarantee that this would also work in a Windows Subsystem for Linux (WSL) .

* `GridapPETSc` default sparse matrix format is 0-based compressed sparse row. This type of sparse matrix storage format can be described by the `SparseMatrixCSR{0,PetscReal,PetscInt}` and `SymSparseMatrixCSR{0,PetscReal,PetscInt}` Julia types as implemented in the [SparseMatricesCSR](https://gridap.github.io/SparseMatricesCSR.jl/stable/) Julia package.
* **When running in MPI parallel mode** (i.e., with a MPI communicator different from `MPI.COMM_SELF`), `GridapPETSc` implements a sort of limited garbage collector in order to automatically deallocate PETSc objects. This garbage collector can be manually triggered by a call to the function `GridapPETSc.gridap_petsc_gc()`. `GridapPETSc` automatically calls this function inside at different strategic points, and **this will be sufficient for most applications**. However, for some applications, with a very frequent allocation of PETSc objects, it might be needed to call this function from application code. This need will be signaled by PETSc via the following internal message error `PETSC ERROR: No more room in array, limit 256 
recompile src/sys/objects/destroy.c with larger value for MAXREGDESOBJS`
## Interesting links

[How to compile PETSc to work with complex arithmetics](https://abhigupta.io/2021/12/08/installing-petsc-complex.html)