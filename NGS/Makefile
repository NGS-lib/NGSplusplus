GXX?=g++

UTIL_DIR=./utility
IO_DIR=./IO
WRITER_DIR=./IO/Writer
PARSER_DIR=./IO/Parser
SRC_DIR=.
DESTDIR=/usr/lib
INCLUDE_DIR=/usr/include/NGS
SAMPLE_DIR=./sample-src
OBJ_DIR=./obj
PROJECT_NAME=libNGS.a
PROJECT_DIR=../libs
NGS_DIR=.
API_DIR=./api
API_ALGORITHMS_DIR=./api/algorithms
API_INTERNAL_DIR=./api/internal
API_BAM_DIR=./api/internal/bam
API_INDEX_DIR=./api/internal/index
API_IO_DIR=./api/internal/io
API_SAM_DIR=./api/internal/sam
API_UTILS_DIR=./api/internal/utils
SHARED_DIR=./shared
BOOST_DIR=./boost

# Nothing should have to be changed after this.
VPATH=UTIL:IO:IO/Parser:IO/Writer
OBJECTS = $(UTIL_DIR)/utility.o
OBJECTS += $(IO_DIR)/uHeader.o $(IO_DIR)/uToken.o 
OBJECTS += $(WRITER_DIR)/uWriter.o $(WRITER_DIR)/uWriterBase.o $(WRITER_DIR)/uWriterBed.o 
OBJECTS += $(WRITER_DIR)/uWriterBed4.o $(WRITER_DIR)/uWriterBed6.o $(WRITER_DIR)/uWriterCustom.o
OBJECTS += $(WRITER_DIR)/uWriterBed3.o
OBJECTS += $(WRITER_DIR)/uWriterBed5.o
OBJECTS += $(WRITER_DIR)/uWriterSam.o $(WRITER_DIR)/uWriterWig.o $(WRITER_DIR)/uWriterBedGraph.o
OBJECTS += $(WRITER_DIR)/uWriterGFF.o $(WRITER_DIR)/uWriterGTF.o
OBJECTS += $(PARSER_DIR)/uParser.o $(PARSER_DIR)/uParserBase.o $(PARSER_DIR)/uParserWig.o 
OBJECTS += $(PARSER_DIR)/uParserSam.o $(PARSER_DIR)/uParserBed.o $(PARSER_DIR)/uParserCustom.o
OBJECTS += $(PARSER_DIR)/uParserBedGraph.o $(PARSER_DIR)/uParserGTF.o $(PARSER_DIR)/uParserGFF.o
OBJECTS += $(PARSER_DIR)/uParserGenePred.o $(PARSER_DIR)/uParserBAM.o
OBJECTS += $(SRC_DIR)/uBasicNGS.o $(SRC_DIR)/uRegion.o $(SRC_DIR)/uTags.o
OBJECTS += $(SRC_DIR)/uGene.o

# Third-Party
OBJECTS += $(API_DIR)/BamAlignment.o
OBJECTS += $(API_DIR)/BamMultiReader.o
OBJECTS += $(API_DIR)/BamReader.o
OBJECTS += $(API_DIR)/BamWriter.o
OBJECTS += $(API_BAM_DIR)/BamHeader_p.o
OBJECTS += $(API_BAM_DIR)/BamMultiReader_p.o
OBJECTS += $(API_BAM_DIR)/BamRandomAccessController_p.o
OBJECTS += $(API_BAM_DIR)/BamReader_p.o
OBJECTS += $(API_BAM_DIR)/BamWriter_p.o
OBJECTS += $(API_INDEX_DIR)/BamIndexFactory_p.o
OBJECTS += $(API_INDEX_DIR)/BamStandardIndex_p.o
OBJECTS += $(API_INDEX_DIR)/BamToolsIndex_p.o
OBJECTS += $(API_IO_DIR)/BamDeviceFactory_p.o
OBJECTS += $(API_IO_DIR)/BamFile_p.o
OBJECTS += $(API_IO_DIR)/BamFtp_p.o
OBJECTS += $(API_IO_DIR)/BamHttp_p.o
OBJECTS += $(API_IO_DIR)/BamPipe_p.o
OBJECTS += $(API_IO_DIR)/BgzfStream_p.o
OBJECTS += $(API_IO_DIR)/ByteArray_p.o
OBJECTS += $(API_IO_DIR)/HostAddress_p.o
OBJECTS += $(API_IO_DIR)/HostInfo_p.o
OBJECTS += $(API_IO_DIR)/HttpHeader_p.o
OBJECTS += $(API_IO_DIR)/ILocalIODevice_p.o
OBJECTS += $(API_IO_DIR)/RollingBuffer_p.o
OBJECTS += $(API_IO_DIR)/TcpSocket_p.o
OBJECTS += $(API_IO_DIR)/TcpSocketEngine_p.o
OBJECTS += $(API_IO_DIR)/TcpSocketEngine_unix_p.o
OBJECTS += $(API_SAM_DIR)/SamFormatParser_p.o
OBJECTS += $(API_SAM_DIR)/SamFormatPrinter_p.o
OBJECTS += $(API_SAM_DIR)/SamHeaderValidator_p.o
OBJECTS += $(API_UTILS_DIR)/BamException_p.o
OBJECTS += $(API_DIR)/SamHeader.o
OBJECTS += $(API_DIR)/SamProgram.o
OBJECTS += $(API_DIR)/SamProgramChain.o
OBJECTS += $(API_DIR)/SamReadGroup.o
OBJECTS += $(API_DIR)/SamReadGroupDictionary.o
OBJECTS += $(API_DIR)/SamSequence.o
OBJECTS += $(API_DIR)/SamSequenceDictionary.o

OBJPROJ = $(addprefix $(OBJ_DIR)/, $(OBJECTS))
PROJECT = $(PROJECT_DIR)/$(PROJECT_NAME)

CFLAGS  = -O2 -std=c++11 -lz
INCLUDES = -I$(NGS_DIR) 

$(shell mkdir -p $(OBJ_DIR)/$(WRITER_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(PARSER_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(UTIL_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_ALGORITHMS_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_INTERNAL_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_BAM_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_SAM_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_IO_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_UTILS_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_INDEX_DIR))
$(shell mkdir -p $(OBJ_DIR)/$(API_SHARED_DIR))


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
	mkdir -p $(INCLUDE_DIR)/$(API_DIR)
	mkdir -p $(INCLUDE_DIR)/$(API_ALGORITHMS_DIR)
	mkdir -p $(INCLUDE_DIR)/$(SHARED_DIR)
	mkdir -p $(INCLUDE_DIR)/$(API_INDEX_DIR)
	mkdir -p $(INCLUDE_DIR)/$(API_IO_DIR)
	mkdir -p $(INCLUDE_DIR)/$(API_UTILS_DIR)
	mkdir -p $(INCLUDE_DIR)/$(API_SAM_DIR)
	mkdir -p $(INCLUDE_DIR)/$(API_BAM_DIR)
	cp $(SRC_DIR)/*.h $(INCLUDE_DIR)
	cp $(SRC_DIR)/../*.h $(INCLUDE_DIR)/..
	cp $(IO_DIR)/*.h $(INCLUDE_DIR)/$(IO_DIR)
	cp $(WRITER_DIR)/*.h $(INCLUDE_DIR)/$(WRITER_DIR)
	cp $(PARSER_DIR)/*.h $(INCLUDE_DIR)/$(PARSER_DIR)
	cp $(UTIL_DIR)/*.h $(INCLUDE_DIR)/$(UTIL_DIR)
	cp $(API_DIR)/*.h $(INCLUDE_DIR)/$(API_DIR)
	cp $(API_ALGORITHMS_DIR)/*.h $(INCLUDE_DIR)/$(API_ALGORITHMS_DIR)
	cp $(API_INDEX_DIR)/*.h $(INCLUDE_DIR)/$(API_INDEX_DIR)
	cp $(API_IO_DIR)/*.h $(INCLUDE_DIR)/$(API_IO_DIR)
	cp $(API_UTILS_DIR)/*.h $(INCLUDE_DIR)/$(API_UTILS_DIR)
	cp $(API_BAM_DIR)/*.h $(INCLUDE_DIR)/$(API_BAM_DIR)
	cp $(API_SAM_DIR)/*.h $(INCLUDE_DIR)/$(API_SAM_DIR)
	cp $(SHARED_DIR)/*.h $(INCLUDE_DIR)/$(SHARED_DIR)
	cp -r $(BOOST_DIR) $(INCLUDE_DIR)

uninstall:
	rm -f $(DESTDIR)/$(PROJECT_NAME)
	rm -f $(INCLUDE_DIR)/../NGS++.h
	rm -rf $(INCLUDE_DIR)

.PHONY : clean
clean: 
	rm -f $(OBJPROJ) $(PROJECT)
