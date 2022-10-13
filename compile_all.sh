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

d_value=False
while getopts "d" opt; do
  case $opt in
    d)
      d_value=True
      ;;
    \?)
      echo "Invalid option: -$OPTARG. Compiling anyway."
      ;;
  esac
done

cmake -DBOUNDS_CHECKS=$d_value ..
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

cmake -DBOUNDS_CHECKS=$d_value ..
make

cd ..
cd ..

# interpolator build

if [ ! -d build ]
then
  mkdir build
fi

cd build

cmake -DBOUNDS_CHECKS=$d_value ..
make

cd ..





