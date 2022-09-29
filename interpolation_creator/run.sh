#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

real2game_root_dir=/home/max/code/real2GAME # the home directory of real2GAME
model_home_dir=/home/max/code/GAME # the home directory of GAME or L-GAME
model_target_id=1 # the model to which to interpolate; 1: GAME 2: L-GAME
model_source_id=1 # the model to interpolate; 1: ICON-global, 2: GAME, 3: ICON-D2
res_id=5 # resolution ID of the GAME grid
n_layers=26 # number of layers of the GAME grid file to use
ny=35 # number of gridpoints in y-direction of L-GAME
nx=35 # number of gridpoints in x-direction of L-GAME
interpol_exp=2.001 # interpolation exponent
lgame_grid="grid.nc" # grid filename of L-GAME
oro_id=1 # orography ID (relevant only for GAME)
export OMP_NUM_THREADS=4 # number of threads used for OMP parallelization

source .sh/root_script.sh
