#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

real2game_root_dir=/home/max/code/real2GAME # the home directory of real2GAME
model_home_dir=/home/max/code/GAME # the home directory of GAME or L-GAME
model_source_id=0 # 0: global ICON, 1: ICON-D2
model_target_id=0 # 0: GAME 1: L-GAME
nlat=25 # number of latitude points of the L-GAME grid
nlon=25 # number of longitude points of the L-GAME grid
interpol_exp=2.001 # interpolation exponent
oro_id=1 # orography ID
export OMP_NUM_THREADS=6 # number of threads used for OMP parallelization

source .sh/root_script.sh
