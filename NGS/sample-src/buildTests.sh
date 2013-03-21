#!/bin/bash

filename=$1

mkdir -p bin

if [ -e "$filename" ]
then
	echo Building $filename in bin/$(basename ${filename%.*})...
	./Scripts/build.sh $filename
	echo Done!
else
	for dir in $(ls | grep -v bin | grep -v README | grep -v buildTests | grep -v Data | grep -v Makefile)
	do
		cpp_file=$(find $dir/* | grep cpp | head -n1)
		if [ -e "$cpp_file" ]
		then
			echo Building $cpp_file in bin/$(basename ${cpp_file%.*})...
			./Scripts/build.sh $cpp_file
			echo Done!
		fi
	done
fi
