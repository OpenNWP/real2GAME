#!/bin/bash

# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

analysis_year=${BASH_ARGV[4]}
analysis_month=${BASH_ARGV[3]}
analysis_day=${BASH_ARGV[2]}
analysis_hour=${BASH_ARGV[1]}
model_home_directory=${BASH_ARGV[0]}
orography_id=${BASH_ARGV[5]}
background_file=${BASH_ARGV[6]}

bin/formatter $analysis_year $analysis_month $analysis_day $analysis_hour
if [ $? -ne 0 ]
then
echo -e ${RED}Formatter failed.$NC
else
echo "Observations formatted for the data assimilation successfully."
fi
