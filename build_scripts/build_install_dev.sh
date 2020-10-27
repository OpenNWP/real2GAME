aim_dir=~/compiled/ndvar
if [ -d $aim_dir/bin ]
then
rm -r $aim_dir/bin
fi
if [ -d $aim_dir/input ]
then
rm -r $aim_dir/input
fi
if [ -d $aim_dir/run_configs ]
then
rm -r $aim_dir/run_configs
fi
if [ -d build ]
then
rm -r build
fi
mkdir $aim_dir/input
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$aim_dir ..
make
ctest
make install
cp ../run.sh $aim_dir
cp -r ../run_configs $aim_dir
chmod +x $aim_dir/run.sh
cd ..
