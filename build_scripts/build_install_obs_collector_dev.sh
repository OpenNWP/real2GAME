# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

aim_dir=~/compiled/ndvar

if [ ! -d $aim_dir ]
then
echo "$aim_dir does not exist. Aborting. Run build_install_dev.sh first."
exit
fi

if [ ! -d $aim_dir/bin ]
then
echo "$aim_dir/bin does not exist. Aborting. Run build_install_dev.sh first."
exit
fi

if [ -d ../build ]
then
rm -r ../build
fi
mkdir ../build && cd ../build
cmake -DCMAKE_INSTALL_PREFIX=$aim_dir ../formatter
make
ctest
make install
cp ../formatter/run_formatter.sh $aim_dir/bin
cd ..
cp -r obs_collector $aim_dir
