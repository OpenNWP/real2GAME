#!/bin/bash

# This downloads data from the global ICON model.

url_base="https://opendata.dwd.de/weather/nwp/icon/grib"

levels_vector=(1 10 19 27 35 43 51 59 67 75 83 90)

# grid properties
# surface height
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_HSURF.grib2.bz2"
url=$url_base"/"$analysis_hour"/hsurf/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename

# model level heights
for level in ${levels_vector[@]}
do
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_"$level"_HHL.grib2.bz2"
url=$url_base"/"$analysis_hour"/hhl/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename
done

# meteorological fields

# surface pressure
filename="icon_global_icosahedral_single-level_"$analysis_year$analysis_month$analysis_day$analysis_hour"_000_PS.grib2.bz2"
url=$url_base"/"$analysis_hour"/ps/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename

# loop over desired model levels
for level in ${levels_vector[@]}
do
# temperature on model levels
filename="icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hour"_000_"$level"_T.grib2.bz2"
url=$url_base"/"$analysis_hour"/t/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename
# u wind on model levels
filename="icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hour"_000_"$level"_U.grib2.bz2"
url=$url_base"/"$analysis_hour"/u/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename
# v wind on model levels
filename="icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hour"_000_"$level"_V.grib2.bz2"
url=$url_base"/"$analysis_hour"/v/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename
# specific humidity on pressure levels
filename="icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hour"_000_"$level"_QV.grib2.bz2"
url=$url_base"/"$analysis_hour"/qv/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename
done









