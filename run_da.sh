#!/bin/bash

# This source file is part of GAME-DA, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME-DA

analysis_year=${BASH_ARGV[4]}
analysis_month=${BASH_ARGV[3]}
analysis_day=${BASH_ARGV[2]}
analysis_hour=${BASH_ARGV[1]}
model_home_directory=${BASH_ARGV[0]}
orography_id=${BASH_ARGV[5]}
background_file=${BASH_ARGV[6]}
game_da_root_dir=${BASH_ARGV[7]}
orography_layers=${BASH_ARGV[8]}
toa=${BASH_ARGV[9]}
omp_num_threads=${BASH_ARGV[10]}

# parallelization
export OMP_NUM_THREADS=$omp_num_threads # relevant only for OMP

echo "This is GAME-DA."
echo "Starting the assimilation process ..."
echo "analysis year: "$analysis_year
echo "analysis month: "$analysis_month
echo "analysis day: "$analysis_day
echo "analysis hour: "$analysis_hour
$game_da_root_dir/build/game-da $analysis_year $analysis_month $analysis_day $analysis_hour $model_home_directory $orography_id $background_file $game_da_root_dir $orography_layers $toa
if [ $? -ne 0 ]
then
echo -e ${RED}Data assimilation failed.$NC
else
echo "Model input file created sucessfully."
fi
