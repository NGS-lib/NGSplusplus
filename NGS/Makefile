GXX?=g++

UTIL_DIR=./utility
IO_DIR=./IO
WRITER_DIR=./IO/Writer
PARSER_DIR=./IO/Parser
SRC_DIR=.
BOOST_DIR=./boost-include
DESTDIR=/usr/lib
INCLUDE_DIR=/usr/include/NGS
SAMPLE_DIR=./sample-src
OBJ_DIR=./obj
PROJECT_NAME=libgeneClass.a
PROJECT_DIR=../libs

# Nothing should have to be changed after this.
VPATH=UTIL:IO:IO/Parser:IO/Writer
OBJECTS = $(UTIL_DIR)/utility.o
OBJECTS += $(IO_DIR)/uHeader.o $(IO_DIR)/uToken.o 
OBJECTS += $(WRITER_DIR)/uWriter.o $(WRITER_DIR)/uWriterBase.o $(WRITER_DIR)/uWriterBed.o 
OBJECTS += $(WRITER_DIR)/uWriterBed4.o $(WRITER_DIR)/uWriterBed6.o $(WRITER_DIR)/uWriterCustom.o
OBJECTS += $(PARSER_DIR)/uParser.o $(PARSER_DIR)/uParserBase.o $(PARSER_DIR)/uParserWig.o 
OBJECTS += $(PARSER_DIR)/uParserSam.o $(PARSER_DIR)/uParserBed.o $(PARSER_DIR)/uParserCustom.o
OBJECTS += $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uGeneException.o $(SRC_DIR)/uRegion.o $(SRC_DIR)/uTags.o
OBJPROJ = $(addprefix $(OBJ_DIR)/, $(OBJECTS))
PROJECT = $(PROJECT_DIR)/$(PROJECT_NAME)

CFLAGS  = -O3 -std=c++11 
INCLUDES = -I$(UTIL_DIR) -I$(BOOST_DIR)

$(shell mkdir -p $(OBJ_DIR)/$(WRITER_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(PARSER_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(UTIL_DIR))

$(PROJECT): $(OBJPROJ)
	mkdir -p $(PROJECT_DIR)
	ar -r -s $(PROJECT) $(OBJPROJ)

$(OBJ_DIR)/%.o: %.cpp
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c -o $@ $<

install:
	cp $(PROJECT) $(DESTDIR)
	mkdir -p $(INCLUDE_DIR)/$(WRITER_DIR)
	mkdir -p $(INCLUDE_DIR)/$(PARSER_DIR)
	mkdir -p $(INCLUDE_DIR)/$(UTIL_DIR)
	cp $(SRC_DIR)/*.h $(INCLUDE_DIR)
	cp $(IO_DIR)/*.h $(INCLUDE_DIR)/$(IO_DIR)
	cp $(WRITER_DIR)/*.h $(INCLUDE_DIR)/$(WRITER_DIR)
	cp $(PARSER_DIR)/*.h $(INCLUDE_DIR)/$(PARSER_DIR)
	cp $(UTIL_DIR)/*.h $(INCLUDE_DIR)/$(UTIL_DIR)

uninstall:
	rm $(DESTDIR)/$(PROJECT_NAME)
	rm -rf $(INCLUDE_DIR)

samples: 
	mkdir -p $(SAMPLE_DIR)/bin
	$(GXX) -Wall -fexceptions -std=c++11 -pg -g -I $(SRC_DIR) -I $(SRC_DIR)/boost-include -c \
	$(SAMPLE_DIR)/DTools/DTools.cpp -o $(SAMPLE_DIR)/DTools/DTools.o 
	$(GXX) -o $(SAMPLE_DIR)/bin/DTools $(SAMPLE_DIR)/DTools/DTools.o $(PROJECT)

.PHONY : clean
clean: 
	rm $(OBJPROJ) $(PROJECT)
