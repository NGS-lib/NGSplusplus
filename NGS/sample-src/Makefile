GXX?=g++
FLAGS=-O2 -Wall -std=c++11
ZLIB=-lz
INCLUDE=-I../../NGS -I../../
LIB=../../libs/libNGS.a

TARGETS = bedToBedgraph/bedTobedgraph
TARGETS += ConvertTool/ConvertTool
TARGETS += CustomFormats/CustomFormats
TARGETS += CutDuplicates/CutDuplicates
TARGETS += DTools/DTools
TARGETS += macsGenePredQuery/MacsGenePredQuery
TARGETS += recurseBed/recurseBedOverlap
TARGETS += ReturnOverlapping/returnOverlapping
TARGETS += splitToChrom/splitToChrom
TARGETS += tutorial/firstProgram
TARGETS += tutorial/linkSample

all: $(TARGETS)

bedToBedgraph/bedTobedgraph: bedToBedgraph/bedTobedgraph.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

ConvertTool/ConvertTool: ConvertTool/ConvertTool.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

CustomFormats/CustomFormats: CustomFormats/CustomFormats.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

CutDuplicates/CutDuplicates: CutDuplicates/CutDuplicates.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

DTools/DTools: DTools/DTools.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

macsGenePredQuery/MacsGenePredQuery: macsGenePredQuery/MacsGenePredQuery.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

recurseBed/recurseBedOverlap: recurseBed/recurseBedOverlap.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

ReturnOverlapping/returnOverlapping: ReturnOverlapping/returnOverlapping.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

splitToChrom/splitToChrom: splitToChrom/splitToChrom.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

tutorial/firstProgram: tutorial/firstProgram.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

tutorial/linkSample: tutorial/linkSample.cpp
	$(GXX) $(FLAGS) $< -o $@ $(LIB) $(INCLUDE) $(ZLIB)

.PHONY: clean
clean:
	rm -f $(TARGETS)
