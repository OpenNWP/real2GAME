#!/bin/bash

if [ -d build ]
then
rm -r build/*
fi

if [ -d formatter/build ]
then
rm -r formatter/build/*
fi
