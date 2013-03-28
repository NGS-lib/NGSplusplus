#!/bin/bash

usage(){
	echo "Usage to build all the samples:"
	echo "	./buildTests.sh"
	echo ""
	echo "Usage to built a specific sample:"
	echo "	./buildTests.sh <filename.cpp>"
	echo ""
	echo "Note: Must be called from the sample-src directory."
	echo ""
}

if [ "$#" -gt 1 ]
then
	usage
elif [ "$#" == 1 ]
then
	filename=$1
	if [ -e "$filename" ]
	then
		echo Building $filename
		./Scripts/build.sh $filename
		echo Done!
	else
		echo "Invalid file name: $filename"
		usage
	fi
else
	for file in $(find * | grep '.cpp')
	do
		echo Building $file
		./Scripts/build.sh $file
	done
fi
