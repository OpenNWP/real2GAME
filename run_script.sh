#!/bin/bash

# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

year=${BASH_ARGV[0]}
month=${BASH_ARGV[1]}
day=${BASH_ARGV[2]}
hour=${BASH_ARGV[3]}
model_home_directory=${BASH_ARGV[4]}

cycle=(0 6 12 18) # the UTC times of the analyses

model_home_directory=/home/max/compiled/game_dev # The directory of the model.

# That's it, now the assimilation process will be started.
source bin/run.sh
