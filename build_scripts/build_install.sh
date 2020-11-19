# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

aim_dir=~/compiled/ndvar
if [ -d $aim_dir ]
then
rm -r $aim_dir
fi
mkdir $aim_dir
if [ -d $aim_dir/bin ]
then
rm -r $aim_dir/bin
fi
if [ -d $aim_dir/input ]
then
rm -r $aim_dir/input
fi
if [ -d ../build ]
then
rm -r ../build
fi
mkdir $aim_dir/input
mkdir ../build && cd ../build
cmake -DCMAKE_INSTALL_PREFIX=$aim_dir ../core
make
ctest
make install
cp ../core/run_ndvar.sh $aim_dir/bin
cd ..
