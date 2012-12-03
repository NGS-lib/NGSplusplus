GXX?=g++

UTIL_DIR=.././utility
IO_DIR=./IO
WRITER_DIR=./IO/Writer
PARSER_DIR=./IO/Parser
SRC_DIR=.
BOOST_DIR=.././boost-include
DESTDIR=./bin
INCLUDE_DIR=../../
SAMPLE_DIR=./sample-src
OBJ_DIR=./obj
LIB_DIR=../../libs
TARGET_DTOOLS=DensityTool
TARGET_CTOOLS=ConversionTool


# Nothing should have to be changed after this.
CFLAGS  = -O3 -std=c++11 
INCLUDES = -I$(UTIL_DIR) -I$(BOOST_DIR)

$(shell mkdir -p $(DESTDIR))
$(shell mkdir -p $(OBJ_DIR))

all: ./DTools/DTools.cpp
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c -o $(OBJ_DIR)/DTools.o ./DTools/DTools.cpp
	$(GXX) $(OBJ_DIR)/DTools.o $(LIB_DIR)/libNGS.a -o $(DESTDIR)/$(TARGET_DTOOLS)
	$(GXX) -c $(CFLAGS) $(INCLUDES) -c -o $(OBJ_DIR)/CTools.o ./ConvertTool/ConvertTool.cpp
	$(GXX) $(OBJ_DIR)/CTools.o $(LIB_DIR)/libNGS.a -o $(DESTDIR)/$(TARGET_CTOOLS)

.PHONY : clean
clean: 
	rm $(OBJPROJ) $(PROJECT)