#!/bin/bash

rm -rf build
mkdir -p build
cd build || exit 1

# recompile LAMMPS without MPI
# cmake -C ../cmake/presets/basic.cmake -D PKG_MBX=yes -D PKG_EXTRA-PAIR=yes -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_C_COMPILER=mpicc  ../cmake
cmake -C ../cmake/presets/basic.cmake -D PKG_MBX=yes ../cmake
make -j2 && make install
cd ..
