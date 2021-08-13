# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/ndvar


if [ ! -d build ]
then
mkdir build
fi

cd build

cmake ..
make

cd ..
