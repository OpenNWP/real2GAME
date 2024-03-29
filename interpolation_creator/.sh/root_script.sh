#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

# sanity checks
if [ $model_source_id -le 0 ] || [ $model_source_id -ge 4 ]
then
  echo "It must be 1 <= model_src_id <= 3."
  echo "Aborting."
  exit
fi
if [ $model_target_id -le 0 ] || [ $model_target_id -ge 3 ]
then
  echo "It must be 1 <= model_target_id <= 2."
  echo "Aborting."
  exit
fi
if [ $model_target_id -eq 1 ] && [ $model_source_id -ne 1 ]
then
  echo "GAME can only be initialized based on ICON-Global."
  echo "Aborting."
  exit
fi

# Firstly, we have to download the grid data.

# global ICON grid data download
if [ $model_source_id -eq 1 ]
then

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
fi

# ICON-D2 grid data download
if [ $model_source_id -eq 3 ]
then

  # properties of the ICON-D2 model
  analysis_delay_min=70
  cycle=(0 3 6 9 12 15 18 21) # the UTC times of the analyses
  url_base="https://opendata.dwd.de/weather/nwp/icon-d2/grib"
  
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
  filename="icon-d2_germany_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_000_0_clat.grib2.bz2"
  url=$url_base"/"$analysis_hour"/clat/"$filename
  wget -q --directory-prefix=$real2game_root_dir/interpolation_creator $url
  bzip2 -d $real2game_root_dir/interpolation_creator/$filename
  filename="icon-d2_germany_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_000_0_clon.grib2.bz2"
  url=$url_base"/"$analysis_hour"/clon/"$filename
  wget -q --directory-prefix=$real2game_root_dir/interpolation_creator $url
  bzip2 -d $real2game_root_dir/interpolation_creator/$filename

  echo "ICON grid data downloaded."
fi

# downloading the SST grid data
echo "Downloading SST grid data ..."
now=$(date +%s)
save_avail_time=$(bc <<< "$now - 1440*60")

ana_day=$(date --utc -d @$save_avail_time +%Y%m%d)

basic_url=https://nomads.ncep.noaa.gov/pub/data/nccf/com/nsst/prod/nsst.
filename=rtgssthr_grb_0.5.grib2
url=$basic_url$ana_day
url=$url/$filename
wget -q --directory-prefix=$real2game_root_dir/interpolation_creator $url
echo "SST grid data downloaded."

# creating the directory for the interpolation files if it does not exist
if [ ! -d $real2game_root_dir/interpolation_files ]
then
  mkdir $real2game_root_dir/interpolation_files
fi

# Now we can execute the interpolation creator itself.
./build/interpolation_creator $analysis_year $analysis_month $analysis_day $analysis_hour $real2game_root_dir $model_home_dir $oro_id $model_target_id $ny $nx $interpol_exp $lgame_grid $res_id $n_layers $model_source_id

# deleting the ICON grid data
rm icon*
rm rtgssthr_grb_0.5.grib2*








