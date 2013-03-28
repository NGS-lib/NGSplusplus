#!/bin/bash

GXX=$1
filename=$2
output=${filename%.*}

$GXX -O2 -Wall $filename -o $output -std=c++11 ../../libs/libNGS.a -I../../NGS -lz
