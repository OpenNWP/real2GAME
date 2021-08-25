#!/bin/bash

# This source file is part of GAME-DA, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME-DA

analysis_year=${BASH_ARGV[3]}
analysis_month=${BASH_ARGV[2]}
analysis_day=${BASH_ARGV[1]}
analysis_hour=${BASH_ARGV[0]}
game_da_root_dir=${BASH_ARGV[4]}

$game_da_root_dir/formatter/build/formatter $analysis_year $analysis_month $analysis_day $analysis_hour $game_da_root_dir
if [ $? -ne 0 ]
then
echo -e ${RED}Formatter failed.$NC
else
echo "Observations formatted for the data assimilation successfully."
fi
