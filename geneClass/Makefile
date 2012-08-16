UTIL_DIR=./utility
#SRC_DIR=./geneClass
SRC_DIR=.
BOOST_DIR=./boost-include

#ar rcs my_library.a file1.o file2.o
###
PROJECT = ../libs/libgeneClass.a
OBJECTS = $(SRC_DIR)/uTags.o $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uFunctions.o $(SRC_DIR)/uRegion.o
OBJECTS += $(SRC_DIR)/utility.o $(SRC_DIR)/uToken.o
UTILS = $(UTIL_DIR)/utility.o
CFLAGS  = -O2 -std=c++11
INCLUDES = -I$(UTIL_DIR) -I$(BOOST_DIR)

all: $(PROJECT) $(UTILS)

.cpp.o:
	/Users/charles/Documents/gcc/dist/bin/g++ -c $(CFLAGS) $(INCLUDES) $<

$(PROJECT): $(OBJECTS)
	libtool -o $(PROJECT) -static $(OBJECTS) $(SRC_DIR)/utility.o
