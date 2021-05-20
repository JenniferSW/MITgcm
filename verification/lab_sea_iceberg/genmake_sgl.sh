#! /usr/bin/env bash

# Written:     2013-10-15
# Last change: 2014-02-11

experiment="lab_sea_iceberg"

cd build

~/MITgcm/tools/genmake2 -rootdir=/home/jwesselbaum/MITgcm \
-mods=/home/jwesselbaum/MITgcm/verification/$experiment/code \
-devel \
-of \
/home/jwesselbaum/MITgcm/tools/build_options/linux_amd64_gfortran_path

make depend
make
