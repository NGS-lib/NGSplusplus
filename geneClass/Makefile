UTIL_DIR=./utility
IO_DIR=./IO
SRC_DIR=.
BOOST_DIR=./boost-include

PROJECT = ../libs/libgeneClass.a
OBJECTS = $(SRC_DIR)/uTags.o $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uFunctions.o $(SRC_DIR)/uRegion.o
OBJECTS += $(SRC_DIR)/uToken.o $(SRC_DIR)/uParser.o $(SRC_DIR)/uHeader.o
UTILS = $(UTIL_DIR)/utility.o
CFLAGS  = -O2 -std=c++11
INCLUDES = -I$(UTIL_DIR) -I$(BOOST_DIR)

all: $(UTILS) $(PROJECT)

.cpp.o:
	g++-4.7 -c $(CFLAGS) $(INCLUDES) $<

$(PROJECT): $(OBJECTS)
	libtool -o $(PROJECT) -static $(OBJECTS) $(SRC_DIR)/utility.o
