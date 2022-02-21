#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

analysis_year=${BASH_ARGV[4]}
analysis_month=${BASH_ARGV[3]}
analysis_day=${BASH_ARGV[2]}
analysis_hour=${BASH_ARGV[1]}
model_home_directory=${BASH_ARGV[0]}
orography_id=${BASH_ARGV[5]}
background_file=${BASH_ARGV[6]}
real2game_home_dir=${BASH_ARGV[7]}
omp_num_threads=${BASH_ARGV[8]}

# parallelization
export OMP_NUM_THREADS=$omp_num_threads # relevant only for OMP

echo "This is real2GAME."

# executing the downloader ...
echo "Starting to download initial data ..."
$real2game_home_dir/downloader/run.sh $real2game_home_dir $analysis_year $analysis_month $analysis_day $analysis_hour
echo "Collection of initial data completed."

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
$real2game_home_dir/build/real2game $analysis_year $analysis_month $analysis_day $analysis_hour $model_home_directory $orography_id $background_file $real2game_home_dir
if [ $? -ne 0 ]
then
echo -e ${RED}real2GAME failed.$NC
else
echo "Model input file created sucessfully."
fi




