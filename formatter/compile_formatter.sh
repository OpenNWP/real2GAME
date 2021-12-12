# This source file is part of GAME-DA, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME-DA


if [ ! -d build ]
then
mkdir build
fi

cd build

cmake ..
make

cd ..
