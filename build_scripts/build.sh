aim_dir=~/compiled/ndvar
if [ -d $aim_dir/input ]
then
rm -r $aim_dir/input
fi
if [ -d build ]
then
rm -r build
fi
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$aim_dir ../ndvar
make
ctest
