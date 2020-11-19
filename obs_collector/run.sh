#!/bin/bash

# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

bin/formatter $analysis_year $analysis_month $analysis_day $analysis_hour
if [ $? -ne 0 ]
then
echo -e ${RED}Formatter failed.$NC
else
echo "Observations formatted for the data assimilation sucessfully."
fi
