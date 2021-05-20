#!/bin/bash

WORKDIR=$(pwd)
builddir=$1
rundir=$2
inputdir=$3


echo "# Locating needed directories"

if [ ! -d $builddir ]; then

  echo "= Directories not found!"
  echo "> Go to parent directory"

  WORKDIR="$WORKDIR/.."
  cd $WORKDIR

  if [ ! -d $builddir ]; then
    exit 1
  fi
fi


echo "= Required directory found."
echo "= Working now in $(pwd)"

if [ ! -d $rundir ]; then
 echo "= run directory not found!"
 echo "ABNORMAL END"
 exit 1
fi

if [ ! -d $inputdir ]; then
 echo "= input directory not found!"
 echo "ABNORMAL END"
 exit 1
fi


#
# Link input and build artifact
#

echo "> Link inputs and artifacts"
cd $rundir
rm *
ln -s ../$inputdir/* .
ln -s ../$builddir/mitgcmuv .

echo "= Links created"
echo "> Run model"


time ./mitgcmuv > output.txt
echo "= Run finished"

