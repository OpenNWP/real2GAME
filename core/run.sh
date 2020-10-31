#!/bin/bash

# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

echo "This is ndvar."
echo "Copyright (C) 2020 The ndvar development team."
echo "Starting the assimilation process ..."
echo "analysis year: "$analysis_year
echo "analysis month: "$analysis_month
echo "analysis day: "$analysis_day
echo "analysis hour: "$analysis_hour
mpirun -np 1 bin/ndvar $analysis_year $analysis_month $analysis_day $analysis_hour $model_home_directory $orography_id
if [ $? -ne 0 ]
then
echo -e ${RED}Data assimilation failed.$NC
else
echo "Model input file created sucessfully."
fi
