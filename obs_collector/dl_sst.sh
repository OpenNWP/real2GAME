#!/bin/bash

# SST downloader

now=$(date +%s)
save_avail_time=$(bc <<< "$now - 1440*60")
ana_day=$(date --utc -d @$save_avail_time +%Y%m%d)

basic_url=https://ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/sst.
filename=rtgssthr_grb_0.083.grib2
url=$basic_url$ana_day
url=$url/$filename
wget --no-verbose --directory-prefix=$game_da_root_dir/input $url
