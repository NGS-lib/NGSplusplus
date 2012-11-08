GXX?=g++

UTIL_DIR=./utility
IO_DIR=./IO
WRITER_DIR=./IO/Writer
PARSER_DIR=./IO/Parser
SRC_DIR=.
BOOST_DIR=./boost-include
DESTDIR=/usr/lib
INCLUDE_DIR=/usr/include
SAMPLE_DIR=./sample-src
OBJ_DIR=./obj
VPATH=UTIL:IO:IO/Parser:IO/Writer

# Nothing should have to be changed after this.
OBJECTS = $(UTIL_DIR)/utility.o
OBJECTS += $(IO_DIR)/uHeader.o $(IO_DIR)/uToken.o 
OBJECTS += $(WRITER_DIR)/uWriter.o $(WRITER_DIR)/uWriterBase.o $(WRITER_DIR)/uWriterBed.o 
OBJECTS += $(WRITER_DIR)/uWriterBed4.o $(WRITER_DIR)/uWriterBed6.o $(WRITER_DIR)/uWriterCustom.o
OBJECTS += $(PARSER_DIR)/uParser.o $(PARSER_DIR)/uParserBase.o $(PARSER_DIR)/uParserWig.o 
OBJECTS += $(PARSER_DIR)/uParserSam.o $(PARSER_DIR)/uParserBed.o $(PARSER_DIR)/uParserCustom.o
OBJECTS += $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uGeneException.o $(SRC_DIR)/uRegion.o $(SRC_DIR)/uTags.o
OBJPROJ = $(addprefix $(OBJ_DIR)/, $(OBJECTS))
PROJECT = ../libs/libgeneClass.a

CFLAGS  = -O3 -std=c++11 
INCLUDES = -I$(UTIL_DIR) -I$(BOOST_DIR)

$(shell mkdir -p $(OBJ_DIR)/$(WRITER_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(PARSER_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(UTIL_DIR))

$(PROJECT): $(OBJPROJ)
	ar -r -s $(PROJECT) $(OBJPROJ)

$(OBJ_DIR)/%.o: %.cpp
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c -o $@ $<

install:
	cp $(PROJECT) $(DESTDIR)
	cp $(SRC_DIR)/*.h $(IO_DIR)/*.h $(WRITER_DIR)/*.h $(PARSER_DIR)/*.h $(INCLUDE_DIR)

samples: 
	mkdir -p $(SAMPLE_DIR)/bin
	$(GXX) -Wall -fexceptions -std=c++11 -pg -g -I $(SRC_DIR) -I $(SRC_DIR)/boost-include -c \
	$(SAMPLE_DIR)/DTools/DTools.cpp -o $(SAMPLE_DIR)/DTools/DTools.o 
	$(GXX) -o $(SAMPLE_DIR)/bin/DTools $(SAMPLE_DIR)/DTools/DTools.o $(PROJECT)

PHONY : clean
clean: 
	rm $(OBJPROJ) $(PROJECT)
