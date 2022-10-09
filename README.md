# fem-tddft

:construction: Work In Progress

**fem-tddft** is a C++ code for materials modeling from first principles using
Kohn-Sham time-dependent density functional theory base on finite-element method
using DEAL II library.

## Build

### Build PETSc

```bash
./configure -PETSC_ARCH=arch-mpi-shared-debug \
            --with-openmpi \
            --with-openblas=1 \
            --with-shared=1 \
            --with-debugging=1
```

### Build SLEPc

Switch to `release` branch

```bash
git checkout release
```
 
```bash
PETSC_DIR=~/software/petsc PETSC_ARCH='arch-mpi-shared-opt' \
./configure --prefix=/opt/slepc/mpi-shared-opt --with-clean --with-precision=single
```

### Build DEAL II

```bash
mkdir build && cd build
```

```bash
cmake .. \
-D DEAL_II_WITH_MPI=ON \
-D DEAL_II_WITH_BOOST=ON -D BOOST_DIR=/usr/include/boost \
-D DEAL_II_WITH_PETSC=ON -D PETSC_DIR=~/software/petsc -D PETSC_ARCH=arch-mpi-shared-opt \
-D DEAL_II_WITH_SLEPC=ON -D SLEPC_DIR=~/software/slepc \
```

```bash
make -j7
```

### Build deal2lkit

When setting up tests, `DEAL_II_DIR` is looked up in the `ENV`, and not in the config!
Tests also required `numdiff` installed.

```bash
cmake .. \
-D DEAL_II_DIR=~/software/dealii-9.2.0/build \
-D D2K_COMPONENT_DOCUMENTATION=OFF \
-D D2K_COMPONENT_EXAMPLES=ON \
-D D2K_ENABLE_TESTING=ON
```

