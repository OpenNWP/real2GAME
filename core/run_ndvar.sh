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

echo "This is ndvar."
echo "Copyright (C) 2020 The ndvar development team."
echo "Starting the assimilation process ..."
echo "analysis year: "$analysis_year
echo "analysis month: "$analysis_month
echo "analysis day: "$analysis_day
echo "analysis hour: "$analysis_hour
mpirun -np 1 bin/ndvar $analysis_year $analysis_month $analysis_day $analysis_hour $model_home_directory $orography_id $background_file
if [ $? -ne 0 ]
then
echo -e ${RED}Data assimilation failed.$NC
else
echo "Model input file created sucessfully."
fi
