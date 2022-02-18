#!/bin/bash

url_base="https://opendata.dwd.de/weather/nwp/icon/grib"

levels_vector=(24 37 50 63 77 90)

# grid properties
# surface height
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_HSURF.grib2.bz2"
url=$url_base"/"$analysis_hour"/hsurf/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename

# horizontal coordinates of the cells
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_CLON.grib2.bz2"
url=$url_base"/"$analysis_hour"/clon/"$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
bzip2 -d $real2game_root_dir/input/$filename
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hour"_CLAT.grib2.bz2"
url=$url_base"/"$analysis_hour"/clat/"$filename
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









