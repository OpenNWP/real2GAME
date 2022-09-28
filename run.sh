#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

analysis_year=${BASH_ARGV[4]}
analysis_month=${BASH_ARGV[3]}
analysis_day=${BASH_ARGV[2]}
analysis_hour=${BASH_ARGV[1]}
model_home_directory=${BASH_ARGV[0]}
oro_id=${BASH_ARGV[5]}
background_file=${BASH_ARGV[6]}
real2game_home_dir=${BASH_ARGV[7]}
omp_num_threads=${BASH_ARGV[8]}
res_id=${BASH_ARGV[9]}
n_layers=${BASH_ARGV[10]}
nsoillays=${BASH_ARGV[11]}
model_src_id=${BASH_ARGV[12]} # 1: ICON-Global, 2: GAME, 3: ICON-D2
model_target_id=${BASH_ARGV[13]} # 1: GAME, 2: L-GAME

# sanity checks
if [ $model_src_id -le 0 ] || [ $model_src_id -ge 4 ] then
  echo "It must be 1 <= model_src_id <= 3."
  echo "Aborting."
  exit
fi
if [ $model_target_id -le 0 ] || [ $model_target_id -ge 3 ] then
  echo "It must be 1 <= model_target_id <= 2."
  echo "Aborting."
  exit
fi
if [ $model_target_id -e 1 ] && [ $model_src_id -ne 1 ] then
  echo "GAME can only be initialized based on ICON-Global."
  echo "Aborting."
  exit
fi

# parallelization
export OMP_NUM_THREADS=$omp_num_threads # relevant only for OMP

echo "This is real2GAME."

# executing the downloader ...
echo "Starting to download initial data ..."
$real2game_home_dir/downloader/run.sh $real2game_home_dir $analysis_year $analysis_month $analysis_day $analysis_hour
echo "Collection of initial data completed."

if [ ! -f $real2game_home_dir/formatter/build/formatter ]
then
  echo "Executable formatter/build/formatter does not exist. Compile first."
  echo "Aborting."
  exit
fi

# reformatting
echo "Reformatting the input data ..."
$real2game_home_dir/formatter/build/formatter $analysis_year $analysis_month $analysis_day $analysis_hour $real2game_home_dir
if [ $? -ne 0 ]
then
  echo -e ${RED}Formatter failed.$NC
else
  echo "Data formatted for the interpolation successfully."
fi

# cleaning the input directory
rm $real2game_home_dir/input/*.grib2

echo "Starting the interpolation process ..."
echo "analysis year: "$analysis_year
echo "analysis month: "$analysis_month
echo "analysis day: "$analysis_day
echo "analysis hour: "$analysis_hour

if [ ! -f $real2game_home_dir/build/real2game ]
then
  echo "Executable /build/real2game does not exist. Compile first."
  echo "Aborting."
  exit
fi

$real2game_home_dir/build/real2game $res_id $n_layers $nsoillays $analysis_year $analysis_month $analysis_day $analysis_hour $model_home_directory $oro_id $background_file $real2game_home_dir
if [ $? -ne 0 ]
then
  echo -e ${RED}real2GAME failed.$NC
else
  echo "Model input file created sucessfully."
fi

# cleaning the input directory
rm $real2game_home_dir/input/*




