#!/usr/bin/env bash

SCRIPT=/home/puck/bin/gif.pml

WORK_DIR=`pwd`
INPUT=$1
OUTPUT=${1%.xyz}.gif
cp $INPUT /tmp/

cd /tmp/
rm -f /tmp/mov*.png
cp $PWD/$INPUT
pymol -c $INPUT $SCRIPT
cd $WORK_DIR

convert -delay 20 -loop 0 /tmp/mov*.png $WORK_DIR/$OUTPUT

rm -f /tmp/mov*.png
