#!/bin/bash

# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

if [ -d interpolation_creator/build ]
then
  rm -r interpolation_creator/build/*
fi

if [ -d formatter/build ]
then
  rm -r formatter/build/*
fi

if [ -d build ]
then
  rm -r build/*
fi

