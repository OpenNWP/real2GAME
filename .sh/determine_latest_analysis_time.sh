#!/bin/bash

# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

# In this file, the time of the analysis to be executed is determined.

now=$(date +%s)
year=$(date --utc -d @$now +%Y)
month=$(date --utc -d @$now +%m)
day=$(date --utc -d @$now +%d)
hour=$(date --utc -d @$now +%H)
# finding the correct analysis hour
hour=${cycle[-1]}
for i in $(seq 0 1 $((${#cycle[@]} - 2)))
do
if [ ${cycle[$i]} -le $hour ] && [ ${cycle[$(($i + 1))]} -gt $hour ]
then
hour=${cycle[$i]}
fi
done
