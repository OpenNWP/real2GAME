#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

real2game_root_dir=${BASH_ARGV[4]}
analysis_year=${BASH_ARGV[3]}
analysis_month=${BASH_ARGV[2]}
analysis_day=${BASH_ARGV[1]}
analysis_hour=${BASH_ARGV[0]}

# downloading weather data
source $real2game_root_dir/downloader/dl_icon-global.sh

# downloading SST
source $real2game_root_dir/downloader/dl_sst.sh
