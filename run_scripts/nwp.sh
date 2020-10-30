#!/bin/bash

# This source file is part of ndvar, which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

cycle=(0 6 12 18) # the UTC times of the analyses

model_home_directory=/home/max/compiled/game_dev # The directory of the model.

# That's it, now the assimilation process will be started.
source ../.sh/determine_latest_analysis_time.sh
source ../bin/run.sh
