# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME


if [ ! -d build ]
then
mkdir build
fi

cd build

cmake ..
make

cd ..
