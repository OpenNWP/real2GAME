#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

# properties of the ICON model
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

./interpolation_creator $analysis_year $analysis_month $analysis_day $analysis_hour $real2game_root_dir $model_home_dir $oro_id
