UTIL_DIR=./utility
IO_DIR=./IO
WRITER_DIR=./IO/Writer
PARSER_DIR=./IO/Parser
SRC_DIR=.
BOOST_DIR=./boost-include
GXX?=g++
DESTDIR=/usr/lib

PROJECT = ../libs/libgeneClass.a
OBJECTS = $(SRC_DIR)/uTags.o $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uFunctions.o $(SRC_DIR)/uRegion.o 
OBJECTS += $(SRC_DIR)/uHeader.o $(SRC_DIR)/uToken.o 
OBJECTS += $(SRC_DIR)/uWriter.o $(SRC_DIR)/uWriterBase.o $(SRC_DIR)/uWriterBed.o 
OBJECTS += $(SRC_DIR)/uWriterBed4.o $(SRC_DIR)/uWriterBed6.o $(SRC_DIR)/uWriterCustom.o
OBJECTS += $(SRC_DIR)/uParser.o $(SRC_DIR)/uParserBase.o $(SRC_DIR)/uParserWig.o 
OBJECTS += $(SRC_DIR)/uParserSam.o $(SRC_DIR)/uParserBed.o $(SRC_DIR)/uParserCustom.o
OBJECTS += $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uGeneException.o $(SRC_DIR)/uRegion.o $(SRC_DIR)/uTags.o
UTILS = $(UTIL_DIR)/utility.o
IO = $(IO_DIR)/uToken.o $(IO_DIR)/uHeader.o 
WRITER = $(WRITER_DIR)/uWriter.o $(WRITER_DIR)/uWriterBase.o $(WRITER_DIR)/uWriterBed.o 
WRITER += $(WRITER_DIR)/uWriterBed4.o $(WRITER_DIR)/uWriterBed6.o $(WRITER_DIR)/uWriterCustom.o
PARSER = $(PARSER_DIR)/uParser.o $(PARSER_DIR)/uParserBase.o $(PARSER_DIR)/uParserWig.o 
PARSER += $(PARSER_DIR)/uParserSam.o $(PARSER_DIR)/uParserBed.o $(PARSER_DIR)/uParserCustom.o
LIBRARY = $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uGeneException.o $(SRC_DIR)/uRegion.o $(SRC_DIR)/uTags.o

CFLAGS  = -O3 -std=c++11 
INCLUDES = -I$(UTIL_DIR) -I$(BOOST_DIR)

all: $(UTILS) $(IO) $(WRITER) $(PARSER) $(PROJECT) $(LIBRARY)

.cpp.o:
	$(GXX) -c $(CFLAGS) $(INCLUDES) $<

$(PROJECT): $(OBJECTS)
	ar -r -s $(PROJECT) $(OBJECTS)

clean: 
	rm *.o $(PROJECT)
