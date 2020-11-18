#!/bin/bash

url_base="https://opendata.dwd.de/weather/nwp/icon/grib"

levels_vector=(70 80 90)

# grid properties
# surface height
wget --directory-prefix=../input $url_base"/"$analysis_hr"/hsurf/icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_HSURF.grib2.bz2"

# horizontal coordinates of the cells
wget --directory-prefix=../input $url_base"/"$analysis_hr"/clon/icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_CLON.grib2.bz2"
wget --directory-prefix=../input $url_base"/"$analysis_hr"/clat/icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_CLAT.grib2.bz2"

# model level heights
for level in ${levels_vector[@]}
do
wget --directory-prefix=../input $url_base"/"$analysis_hr"/hhl/icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_"$level"_HHL.grib2.bz2"
done

# meteorological fields

# surface pressure
wget --directory-prefix=../input $url_base"/"$analysis_hr"/ps/icon_global_icosahedral_single-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_PS.grib2.bz2"

# temperatures on model levels
for level in ${levels_vector[@]}
do
wget --directory-prefix=../input $url_base"/"$analysis_hr"/t/icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_"$level"_T.grib2.bz2"
done

# specific humidity on pressure levels
for level in ${levels_vector[@]}
do
wget --directory-prefix=../input $url_base"/"$analysis_hr"/qv/icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_"$level"_QV.grib2.bz2"
done
