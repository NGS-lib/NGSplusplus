#!/bin/bash

filename=$1
output=${filename%.*}

g++-4.7 -O2 -Wall $filename -o $output -std=c++11 ../../libs/libNGS.a -I../../NGS -lz
