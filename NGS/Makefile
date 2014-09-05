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
OBJ_DIR_STATIC=$(OBJ_DIR)/static
OBJ_DIR_SHARED=$(OBJ_DIR)/shared
PROJECT_NAME_STATIC=libNGS.a
PROJECT_NAME_SHARED=libNGS.so
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
BOOST_IOSTREAM=./libs/iostreams/src
BB_SRC_DIR=./bbfile/src
BB_INCLUDE_DIR=./bbfile/include/

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
OBJECTS += $(PARSER_DIR)/uParserBedGraph.o $(PARSER_DIR)/uParserGTF.o $(PARSER_DIR)/uParserGFF.o $(PARSER_DIR)/uParserBW.o
OBJECTS += $(PARSER_DIR)/uParserGenePred.o $(PARSER_DIR)/uParserBAM.o
OBJECTS += $(SRC_DIR)/uBasicNGS.o $(SRC_DIR)/uRegion.o $(SRC_DIR)/uTags.o
OBJECTS += $(SRC_DIR)/uGene.o

OBJECTS += $(BB_SRC_DIR)/BBFileReader.o $(BB_SRC_DIR)/BBTotalSummaryBlock.o $(BB_SRC_DIR)/BBZoomLevelFormat.o $(BB_SRC_DIR)/BBZoomLevelHeader.o 
OBJECTS += $(BB_SRC_DIR)/BBZoomLevels.o $(BB_SRC_DIR)/BigWigDataBlock.o $(BB_SRC_DIR)/BigWigIterator.o $(BB_SRC_DIR)/BigWigSection.o 
OBJECTS += $(BB_SRC_DIR)/BigWigSectionHeader.o $(BB_SRC_DIR)/BigWigSectionHelper.o $(BB_SRC_DIR)/BPTree.o $(BB_SRC_DIR)/BPTreeChildNode.o 
OBJECTS += $(BB_SRC_DIR)/BPTreeChildNodeItem.o $(BB_SRC_DIR)/BPTreeHeader.o $(BB_SRC_DIR)/BPTreeLeafNode.o $(BB_SRC_DIR)/BPTreeLeafNodeItem.o 
OBJECTS += $(BB_SRC_DIR)/endian_helper.o $(BB_SRC_DIR)/RPChromosomeRegion.o $(BB_SRC_DIR)/RPTree.o $(BB_SRC_DIR)/RPTreeChildNode.o 
OBJECTS += $(BB_SRC_DIR)/RPTreeChildNodeItem.o $(BB_SRC_DIR)/RPTreeHeader.o $(BB_SRC_DIR)/RPTreeLeafNode.o $(BB_SRC_DIR)/RPTreeLeafNodeItem.o 
OBJECTS += $(BB_SRC_DIR)/RPTreeNodeProxy.o $(BB_SRC_DIR)/ZoomDataRecord.o $(BB_SRC_DIR)/ZoomLevelIterator.o 




# Third-Party
OBJECTS += $(API_DIR)/BamAlignment.o
OBJECTS += $(BOOST_IOSTREAM)/zlib.o
OBJECTS += $(BOOST_IOSTREAM)/gzip.o
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





OBJPROJ_STATIC = $(addprefix $(OBJ_DIR_STATIC)/, $(OBJECTS))
OBJPROJ_SHARED = $(addprefix $(OBJ_DIR_SHARED)/, $(OBJECTS))
PROJECT_STATIC = $(PROJECT_DIR)/$(PROJECT_NAME_STATIC)
PROJECT_SHARED = $(PROJECT_DIR)/$(PROJECT_NAME_SHARED)

CFLAGS_STATIC  =-O2 -std=c++11 -lz
CFLAGS_SHARED = $(CFLAGS_STATIC) -fPIC
LFLAGS =-s -static-libgcc -lz -s
INCLUDES =-I$(NGS_DIR) -I$(BB_INCLUDE_DIR) 

$(shell mkdir -p $(OBJ_DIR_STATIC)/$(WRITER_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(PARSER_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(UTIL_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(BOOST_IOSTREAM))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(BOOST_LIB))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_ALGORITHMS_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_INTERNAL_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_BAM_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_SAM_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_IO_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_UTILS_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_INDEX_DIR))
$(shell mkdir -p $(OBJ_DIR_STATIC)/$(API_SHARED_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(WRITER_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(PARSER_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(UTIL_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(BOOST_LIB))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(BOOST_IOSTREAM))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_ALGORITHMS_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_INTERNAL_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_BAM_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_SAM_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_IO_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_UTILS_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_INDEX_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(API_SHARED_DIR))
$(shell mkdir -p $(OBJ_DIR_SHARED)/$(BB_SRC_DIR))

.PHONY: static
static: $(PROJECT_STATIC)
	ar -r -s $(PROJECT_STATIC) $(OBJPROJ_STATIC)

.PHONY: shared
shared: $(PROJECT_SHARED)
	$(GXX) -shared -o $(PROJECT_SHARED) $(OBJPROJ_SHARED) $(LFLAGS)

$(PROJECT_STATIC): $(OBJPROJ_STATIC)
	mkdir -p $(PROJECT_DIR)

$(PROJECT_SHARED): $(OBJPROJ_SHARED)
	mkdir -p $(PROJECT_DIR)

$(OBJ_DIR_STATIC)/%.o: %.cpp
	$(GXX) -c $(CFLAGS_STATIC) $(INCLUDES) -c -o $@ $<

$(OBJ_DIR_SHARED)/%.o: %.cpp
	$(GXX) -c $(CFLAGS_SHARED) $(INCLUDES) -c -o $@ $<

install:
	if [ -f $(PROJECT_STATIC) ]; then cp $(PROJECT_STATIC) $(DESTDIR); fi;		
	if [ -f $(PROJECT_SHARED) ]; then cp $(PROJECT_SHARED) $(DESTDIR); fi;		

	mkdir -p $(INCLUDE_DIR)/$(WRITER_DIR)
	mkdir -p $(INCLUDE_DIR)/$(PARSER_DIR)
	mkdir -p $(INCLUDE_DIR)/$(UTIL_DIR)
	mkdir -p $(INCLUDE_DIR)/$(API_DIR)
	mkdir -p $(INCLUDE_DIR)/$(BOOST_LIB)
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
	rm -f $(DESTDIR)/$(PROJECT_NAME_STATIC)
	rm -f $(DESTDIR)/$(PROJECT_NAME_SHARED)
	rm -f $(INCLUDE_DIR)/../NGS++.h
	rm -rf $(INCLUDE_DIR)

.PHONY : clean
clean: 
	rm -f $(OBJPROJ_STATIC) $(OBJPROJ_SHARED) $(PROJECT_STATIC) $(PROJECT_SHARED)
