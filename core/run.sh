#!/bin/bash

# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

echo "This is ndvar."
echo "Starting the assimilation process ..."
echo "analysis year: "$year
echo "analysis month: "$month
echo "analysis day: "$day
echo "analysis hour: "$hour
mpirun -np 1 ../bin/ndvar $year $month $day $hour $model_home_directory
if [ $? -ne 0 ]
then
echo -e ${RED}Data assimilation failed.$NC
else
echo "Model input file created sucessfully."
fi
