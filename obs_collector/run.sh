#!/bin/bash

analysis_year=${BASH_ARGV[0]}
analysis_month=${BASH_ARGV[1]}
analysis_day=${BASH_ARGV[2]}
analysis_hr=${BASH_ARGV[3]}

# surface height
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/hsurf/icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_HSURF.grib2.bz2"

# surface pressure
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/ps/icon_global_icosahedral_single-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_PS.grib2.bz2"

# model level heights
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/hhl/icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_90_HHL.grib2.bz2"
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/hhl/icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_80_HHL.grib2.bz2"
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/hhl/icon_global_icosahedral_time-invariant_"$analysis_year$analysis_month$analysis_day$analysis_hr"_70_HHL.grib2.bz2"

# temperatures on model levels
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/t/icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_90_T.grib2.bz2"
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/t/icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_80_T.grib2.bz2"
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/t/icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_70_T.grib2.bz2"

# specific humidity on pressure levels
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/qv/icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_90_QV.grib2.bz2"
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/qv/icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_80_QV.grib2.bz2"
wget --directory-prefix=../input "https://opendata.dwd.de/weather/nwp/icon/grib/"$analysis_hr"/qv/icon_global_icosahedral_model-level_"$analysis_year$analysis_month$analysis_day$analysis_hr"_000_70_QV.grib2.bz2"
