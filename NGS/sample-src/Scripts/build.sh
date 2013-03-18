#!/bin/bash

filename=$1
output=bin/$(basename ${filename%.*})

g++-4.7 -lz -O3 -Wall $filename -o $output -std=c++11 ../../libs/libNGS.a -I../../NGS
