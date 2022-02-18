#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

./compile_real2game.sh

cd formatter
./compile_formatter.sh
cd ..
