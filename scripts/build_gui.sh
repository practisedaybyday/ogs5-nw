#!/usr/bin/env bash

# Goto sources directory
cd ..

# Cleanup
rm -rf build_gui

mkdir -p build_gui && cd build_gui
wget http://gist.github.com/raw/527084/c32a0337164f7055ff3b7c92f07f883f4afd3966/cmake.ogs.sh
chmod +x cmake.ogs.sh
./cmake.ogs.sh -DOGS_USE_QT=ON -DOGS_PACKAGING=ON ..
make -j
make package -j
