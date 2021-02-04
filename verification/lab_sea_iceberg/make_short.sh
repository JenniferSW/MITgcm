#!/bin/bash

WORKDIR=$(pwd)

echo "= go to build directory"

if [ ! -d "build" ] ; then
   echo "no build directory found"
   exit 1
fi

cd build

echo "= update compiled files"
make

echo "= leave build directory "
cd ..

