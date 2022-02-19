#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

model_home_directory=/home/max/code/GAME
orography_id=1
background_file=/home/max/code/GAME/standard_oro1.nc
real2game_root_dir=/home/max/code/real2GAME
omp_num_threads=1

# END OF USUAL INPUT SECTION

# setting the time of the analysis
analysis_delay_min=175
cycle=(0 6 12 18) # the UTC times of the analyses
now=$(date +%s)
now=$(($now - $analysis_delay_min*60))
analysis_year=$(date --utc -d @$now +%Y)
analysis_month=$(date --utc -d @$now +%m)
analysis_day=$(date --utc -d @$now +%d)
now_hour=$(date --utc -d @$now +%k)
# finding the correct analysis hour
analysis_hour=${cycle[-1]}
for i in $(seq 0 1 $((${#cycle[@]} - 2)))
do
if [ ${cycle[$i]} -le $now_hour ] && [ ${cycle[$(($i + 1))]} -gt $now_hour ]
then
analysis_hour=${cycle[$i]}
fi
done

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
$real2game_root_dir/downloader/run.sh $real2game_root_dir $analysis_year $analysis_month $analysis_day $analysis_hour_extended_string
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




