#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

d_value=False
f_value=False
while getopts "df" opt; do
  case $opt in
    d) # debugging flag
      d_value=True
      ;;
    f) # aggressive optimization flag
      f_value=True
      ;;
    \?)
      echo "Invalid option: -$OPTARG. Compiling anyway."
      ;;
  esac
done

# interpolation_creator build

cd interpolation_creator
if [ ! -d build ]
then
  mkdir build
fi

cd build

cmake -DDEBUGGING=$d_value ..
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

cmake -DDEBUGGING=$d_value -DFAST=$f_value ..
make

cd ..
cd ..

# interpolator build

if [ ! -d build ]
then
  mkdir build
fi

cd build

cmake -DDEBUGGING=$d_value -DFAST=$f_value ..
make

cd ..











