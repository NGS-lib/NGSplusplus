UTIL_DIR=./utility
IO_DIR=./IO
WRITER_DIR=./IO/Writer
SRC_DIR=.
BOOST_DIR=./boost-include

PROJECT = ../libs/libgeneClass.a
OBJECTS = $(SRC_DIR)/uTags.o $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uFunctions.o $(SRC_DIR)/uRegion.o $(SRC_DIR)/uParser.o
OBJECTS += $(SRC_DIR)/uHeader.o $(SRC_DIR)/uToken.o 
OBJECTS += $(SRC_DIR)/uWriter.o $(SRC_DIR)/uWriterBase.o $(SRC_DIR)/uWriterBed.o
UTILS = $(UTIL_DIR)/utility.o
IO = $(IO_DIR)/uToken.o $(IO_DIR)/uHeader.o 
WRITER = $(WRITER_DIR)/uWriter.o $(WRITER_DIR)/uWriterBase.o $(WRITER_DIR)/uWriterBed.o
CFLAGS  = -O2 -std=c++11 -Wl,--whole-archive -lgeneClass -Wl,--no-whole-archive
INCLUDES = -I$(UTIL_DIR) -I$(BOOST_DIR)

all: $(UTILS) $(IO) $(WRITER) $(PROJECT)

.cpp.o:
	g++-4.7 -c $(CFLAGS) $(INCLUDES) $<

$(PROJECT): $(OBJECTS)
	libtool -o $(PROJECT) -static $(OBJECTS) $(SRC_DIR)/utility.o
