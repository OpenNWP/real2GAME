#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

# Firstly, we have to download the grid data.

# properties of the global ICON model
analysis_delay_min=175
cycle=(0 6 12 18) # the UTC times of the analyses
url_base="https://opendata.dwd.de/weather/nwp/icon/grib"

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

if [ $analysis_hour -lt 10 ]
then
  analysis_hour="0"$analysis_hour
fi

echo "Downloading ICON grid data ..."

# downloading horizontal coordinates of the cells
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_CLAT.grib2.bz2"
url=$url_base"/"$analysis_hour"/clat/"$filename
wget -q --directory-prefix=$real2game_root_dir/interpolation_creator $url
bzip2 -d $real2game_root_dir/interpolation_creator/$filename
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_CLON.grib2.bz2"
url=$url_base"/"$analysis_hour"/clon/"$filename
wget -q --directory-prefix=$real2game_root_dir/interpolation_creator $url
bzip2 -d $real2game_root_dir/interpolation_creator/$filename

echo "ICON grid data downloaded."

# Now we can execute the interpolation creator itself.
./build/interpolation_creator $analysis_year $analysis_month $analysis_day $analysis_hour $real2game_root_dir $model_home_dir $oro_id $model_target_id $nlat $nlon $interpol_exp $lgame_grid

# deleting the ICON grid data
rm icon*










