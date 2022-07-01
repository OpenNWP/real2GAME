#!/bin/bash

# SST downloader

now=$(date +%s)
save_avail_time=$(bc <<< "$now - 1440*60")

ana_day=$(date --utc -d @$save_avail_time +%Y%m%d)

basic_url=https://nomads.ncep.noaa.gov/pub/data/nccf/com/nsst/prod/nsst.
filename=rtgssthr_grb_0.5.grib2
url=$basic_url$ana_day
url=$url/$filename
wget --no-verbose --directory-prefix=$real2game_root_dir/input $url
