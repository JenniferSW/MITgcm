#! /usr/bin/env bash

# Written:     2013-10-15
# Last change: 2014-02-11

experiment="MITgcm/verification/tutorial_barotropic_gyre"

cd build

~/MITgcm_home/MITgcm/tools/genmake2 -rootdir=/home/jwesselbaum/MITgcm_home/MITgcm \
-mods=/home/jwesselbaum/MITgcm_home/$experiment/code \
-of \
/home/jwesselbaum/MITgcm_home/MITgcm/tools/build_options/linux_amd64_gfortran_path

make depend
make