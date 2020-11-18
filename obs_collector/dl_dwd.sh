#!/bin/bash

url_base="https://opendata.dwd.de/weather/nwp/icon/grib"

levels_vector=(70 80 90)

# grid properties
# surface height
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_HSURF.grib2.bz2"
url=$url_base"/"$analysis_hr"/hsurf/"$filename
wget --directory-prefix=../input $url
bzip2 -d ../input/$filename

# horizontal coordinates of the cells
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_CLON.grib2.bz2"
url=$url_base"/"$analysis_hr"/clon/"$filename
wget --directory-prefix=../input $url
bzip2 -d ../input/$filename
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_CLAT.grib2.bz2"
url=$url_base"/"$analysis_hr"/clat/"$filename
wget --directory-prefix=../input $url
bzip2 -d ../input/$filename

# model level heights
for level in ${levels_vector[@]}
do
filename="icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_"$level"_HHL.grib2.bz2"
url=$url_base"/"$analysis_hr"/hhl/"$filename
wget --directory-prefix=../input $url
bzip2 -d ../input/$filename
done

# meteorological fields

# surface pressure
filename="icon_global_icosahedral_single-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_PS.grib2.bz2"
url=$url_base"/"$analysis_hr"/ps/"$filename
wget --directory-prefix=../input $url
bzip2 -d ../input/$filename

# temperatures on model levels
for level in ${levels_vector[@]}
do
filename="icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_"$level"_T.grib2.bz2"
url=$url_base"/"$analysis_hr"/t/"$filename
wget --directory-prefix=../input $url
bzip2 -d ../input/$filename
done

# specific humidity on pressure levels
for level in ${levels_vector[@]}
do
filename="icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_"$level"_QV.grib2.bz2"
url=$url_base"/"$analysis_hr"/qv/"$filename
wget --directory-prefix=../input $url
bzip2 -d ../input/$filename
done




