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
real2game_root_dir=${BASH_ARGV[7]}
omp_num_threads=${BASH_ARGV[8]}

# parallelization
export OMP_NUM_THREADS=$omp_num_threads # relevant only for OMP

echo "This is real2GAME."
echo "Starting the interpolation process ..."
echo "analysis year: "$analysis_year
echo "Analysis month: "$analysis_month
echo "Analysis day: "$analysis_day
echo "Analysis hour: "$analysis_hour
$real2game_root_dir/build/real2game $analysis_year $analysis_month $analysis_day $analysis_hour $model_home_directory $orography_id $background_file $real2game_root_dir
if [ $? -ne 0 ]
then
echo -e ${RED}Data assimilation failed.$NC
else
echo "Model input file created sucessfully."
fi
