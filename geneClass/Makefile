UTIL_DIR=./utility
IO_DIR=./IO
WRITER_DIR=./IO/Writer
PARSER_DIR=./IO/Parser
SRC_DIR=.
BOOST_DIR=./boost-include
GXX?=g++
DESTDIR=/usr/lib

SAMPLE_DIR=./sample-src

PROJECT = ../libs/libgeneClass.a

SOURCE = $(UTIL_DIR)/utility.cpp
SOURCE += $(IO_DIR)/uHeader.cpp $(IO_DIR)/uToken.cpp 
SOURCE += $(WRITER_DIR)/uWriter.cpp $(WRITER_DIR)/uWriterBase.cpp $(WRITER_DIR)/uWriterBed.cpp 
SOURCE += $(WRITER_DIR)/uWriterBed4.cpp $(WRITER_DIR)/uWriterBed6.cpp $(WRITER_DIR)/uWriterCustom.cpp
SOURCE += $(PARSER_DIR)/uParser.cpp $(PARSER_DIR)/uParserBase.cpp $(PARSER_DIR)/uParserWig.cpp 
SOURCE += $(PARSER_DIR)/uParserSam.cpp $(PARSER_DIR)/uParserBed.cpp $(PARSER_DIR)/uParserCustom.cpp
SOURCE += $(SRC_DIR)/uFormatBase.cpp $(SRC_DIR)/uGeneException.cpp $(SRC_DIR)/uRegion.cpp $(SRC_DIR)/uTags.cpp
OBJECTS = $(SOURCE:.cpp=.o)

CFLAGS  = -O3 -std=c++11 
INCLUDES = -I$(UTIL_DIR) -I$(BOOST_DIR)

$(PROJECT): $(OBJECTS)
	ar -r -s $(PROJECT) $(OBJECTS)

%.o: %.cpp
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(UTIL_DIR)/%.o: $(UTIL_DIR)/%.cpp
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(IO_DIR)/%.o: $(IO_DIR)/%.cpp
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(WRITER_DIR)/%.o: $(WRITER_DIR)/%.cpp
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(PARSER_DIR)/%.o: $(PARSER_DIR)/%.cpp
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c o $@ $<

samples: 
	mkdir -p $(SAMPLE_DIR)/bin
	$(GXX) -Wall -fexceptions -std=c++11 -pg -g -I $(SRC_DIR) -I $(SRC_DIR)/boost-include -c \
	$(SAMPLE_DIR)/DTools/DTools.cpp -o $(SAMPLE_DIR)/DTools/DTools.o 
	$(GXX) -o $(SAMPLE_DIR)/bin/DTools $(SAMPLE_DIR)/DTools/DTools.o $(PROJECT)

clean: 
	rm *.o $(IO_DIR)/*.o $(WRITER_DIR)/*.o $(PARSER_DIR)/*.o $(PROJECT)
