#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

analysis_year=2022
analysis_month=02
analysis_day=20
analysis_hour=0
model_home_directory=/home/max/code/GAME
orography_id=1
background_file=/home/max/code/GAME/standard_oro1.nc
real2game_root_dir=/home/max/code/real2GAME
omp_num_threads=2

# analysis hour in a special format
analysis_hour_extended_string=$analysis_hour
if [ $analysis_hour -lt 10 ]
then
analysis_hour_extended_string="0$analysis_hour"
fi

# parallelization
export OMP_NUM_THREADS=$omp_num_threads # relevant only for OMP


echo "This is real2GAME."

# executing the downloader ...
echo "Starting to download initial data ..."
$real2game_home_dir/downloader/run.sh $real2game_home_dir $analysis_year $analysis_month $analysis_day $analysis_hour_extended_string
echo "Collection of initial data completed."

# reformatting
echo "Reformatting the input data .."
$real2game_root_dir/formatter/build/formatter $analysis_year $analysis_month $analysis_day $analysis_hour $real2game_root_dir
if [ $? -ne 0 ]
then
echo -e ${RED}Formatter failed.$NC
else
echo "Data formatted for the interpolation successfully."
fi

echo "Starting the interpolation process ..."
echo "analysis year: "$analysis_year
echo "analysis month: "$analysis_month
echo "analysis day: "$analysis_day
echo "analysis hour: "$analysis_hour
$real2game_root_dir/build/real2game $analysis_year $analysis_month $analysis_day $analysis_hour $model_home_directory $orography_id $background_file $real2game_root_dir
if [ $? -ne 0 ]
then
echo -e ${RED}real2GAME failed.$NC
else
echo "Model input file created sucessfully."
fi




