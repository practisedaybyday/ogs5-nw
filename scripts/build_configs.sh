#!/usr/bin/env bash

# Goto sources directory
cd ..

# Cleanup
rm -rf Release build_fem build_sp build_mpi build_gems build_pqc build_brns build_mkl build_lis

# Executables will be copied to Release/
mkdir -p Release

# Build FEM
mkdir -p build_fem && cd build_fem
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_FEM=ON -DOGS_DONT_USE_QT=ON ..
make -j
cp bin/ogs ../Release/ogs
cd ..

# Build FEM_SP
mkdir -p build_sp && cd build_sp
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_FEM_SP=ON -DOGS_DONT_USE_QT=ON ..
make -j
cp bin/ogs ../Release/ogs_sp
cd ..

# Build FEM_MPI
mkdir -p build_mpi && cd build_mpi
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_FEM_MPI=ON -DOGS_DONT_USE_QT=ON ..
make -j
cp bin/ogs ../Release/ogs_mpi
cd ..

# Build FEM_GEMS
mkdir -p build_gems && cd build_gems
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_FEM_GEMS=ON -DOGS_DONT_USE_QT=ON ..
make -j
cp bin/ogs ../Release/ogs_gems
cd ..

# Build FEM_PQC
mkdir -p build_pqc && cd build_pqc
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_FEM_PQC=ON -DOGS_DONT_USE_QT=ON ..
make -j
cp bin/ogs ../Release/ogs_pqc
cd ..

# Build FEM_BRNS
mkdir -p build_brns && cd build_brns
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_FEM_BRNS=ON -DOGS_DONT_USE_QT=ON ..
make -j
cp bin/ogs ../Release/ogs_brns
cd ..

# Build FEM_MKL
mkdir -p build_mkl && cd build_mkl
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_FEM_MKL=ON -DOGS_DONT_USE_QT=ON ..
make -j
cp bin/ogs ../Release/ogs_mkl
cd ..

# Build FEM_LIS
mkdir -p build_lis && cd build_lis
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_FEM_LIS=ON -DOGS_DONT_USE_QT=ON ..
make -j
cp bin/ogs ../Release/ogs_lis
cd ..
