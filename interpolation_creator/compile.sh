#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

gcc -O2 src/* -fopenmp -lnetcdf -lm -lgeos95 -Wall -o interpolation_creator
