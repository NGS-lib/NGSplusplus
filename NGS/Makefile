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
PROJECT_NAME=libNGS.a
PROJECT_DIR=../libs

# Nothing should have to be changed after this.
VPATH=UTIL:IO:IO/Parser:IO/Writer
OBJECTS = $(UTIL_DIR)/utility.o
OBJECTS += $(IO_DIR)/uHeader.o $(IO_DIR)/uToken.o 
OBJECTS += $(WRITER_DIR)/uWriter.o $(WRITER_DIR)/uWriterBase.o $(WRITER_DIR)/uWriterBed.o 
OBJECTS += $(WRITER_DIR)/uWriterBed4.o $(WRITER_DIR)/uWriterBed6.o $(WRITER_DIR)/uWriterCustom.o
OBJECTS += $(WRITER_DIR)/uWriterSam.o $(WRITER_DIR)/uWriterWig.o $(WRITER_DIR)/uWriterBedGraph.o
OBJECTS += $(WRITER_DIR)/uWriterGFF.o
OBJECTS += $(PARSER_DIR)/uParser.o $(PARSER_DIR)/uParserBase.o $(PARSER_DIR)/uParserWig.o 
OBJECTS += $(PARSER_DIR)/uParserSam.o $(PARSER_DIR)/uParserBed.o $(PARSER_DIR)/uParserCustom.o
OBJECTS += $(PARSER_DIR)/uParserBedGraph.o $(PARSER_DIR)/uParserGTF.o $(PARSER_DIR)/uParserGFF.o
OBJECTS += $(SRC_DIR)/uFormatBase.o $(SRC_DIR)/uRegion.o $(SRC_DIR)/uTags.o
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
	cp $(SRC_DIR)/../*.h $(INCLUDE_DIR)/..
	cp $(IO_DIR)/*.h $(INCLUDE_DIR)/$(IO_DIR)
	cp $(WRITER_DIR)/*.h $(INCLUDE_DIR)/$(WRITER_DIR)
	cp $(PARSER_DIR)/*.h $(INCLUDE_DIR)/$(PARSER_DIR)
	cp $(UTIL_DIR)/*.h $(INCLUDE_DIR)/$(UTIL_DIR)
	cp -r $(BOOST_DIR) $(INCLUDE_DIR)

uninstall:
	rm -f $(DESTDIR)/$(PROJECT_NAME)
	rm -f $(INCLUDE_DIR)/../NGS++.h
	rm -rf $(INCLUDE_DIR)

.PHONY : clean
clean: 
	rm -f $(OBJPROJ) $(PROJECT)
