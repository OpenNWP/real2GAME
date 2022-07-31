#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

# interpolation_creator build

cd interpolation_creator
if [ ! -d build ]
then
  mkdir build
fi

cd build

cmake ..
make

cd ..
cd ..

# formatter build

cd formatter
if [ ! -d build ]
then
  mkdir build
fi

cd build

cmake ..
make

cd ..
cd ..

# interpolator build

if [ ! -d build ]
then
  mkdir build
fi

cd build

cmake ..
make

cd ..
