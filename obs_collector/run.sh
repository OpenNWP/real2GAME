#!/bin/bash

# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/ndvar

analysis_year=${BASH_ARGV[3]}
analysis_month=${BASH_ARGV[2]}
analysis_day=${BASH_ARGV[1]}
analysis_hour=${BASH_ARGV[0]}

source dl_dwd.sh
