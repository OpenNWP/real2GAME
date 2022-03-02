#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

real2game_root_dir=${BASH_ARGV[4]}
analysis_year=${BASH_ARGV[3]}
analysis_month=${BASH_ARGV[2]}
analysis_day=${BASH_ARGV[1]}
analysis_hour=${BASH_ARGV[0]}
model_source_id=${BASH_ARGV[5]}

if [ $model_source_id -eq 0 ]
then
source $real2game_root_dir/downloader/dl_icon-global.sh
fi

if [ $model_source_id -eq 1 ]
then
source $real2game_root_dir/downloader/dl_icon-d2.sh
fi

# downloading SST
source $real2game_root_dir/downloader/dl_sst.sh
