#ifndef BBFILEREADER_H_INCLUDED
#define BBFILEREADER_H_INCLUDED

// ***************************************************************************
//   BBFileReader.h (c) 2014
//   Copyright @ Alexei Nordell-Markovits : Sherbrooke University
//
//    This file is part BWReader library.
//
//    The BWReader library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU  General Public License
//    along with this program (gpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// ***************************************************************************


// The code structure and logic is based on the original IGV BBFileReader.
// The original code can be found here:
//https://github.com/broadinstitute/IGV
// The code was originally released under the LGPL 2.1
//http://www.opensource.org/licenses/lgpl-2.1.php).
// Our thanks to the IGV team for making the code available.


//The following is from the original files of the Broad IGV BBReader
/*
*   Broad Institute Interactive Genome Viewer Big Binary File (BBFile) Reader
*   -   File reader for UCSC BigWig and BigBed file types.
*
*   Notes:   Table entries refer to Jim Kent of UCSC's document description:
*           "BigWig and BigBed: Enabling Browsing of Large Distributed Data Sets",
*           November 2009.
*
*           The overall binary file layout is defined in Table B of the document.
*
*   BBFile Reader sequences through this binary file layout:
*
*   1) Reads in BBFile Header Table C and determine if file is a valid Big Bed or Big Wig
*      binary file type.
*
*   2) Reads in Zoom Header Tables D if zoom data is present, as defined by zoomLevels
*       in Table C, one for each zoom level.
*
*   3) Reads in  the AutoSQL block Table B if present, as referenced by autoSqlOffset in Table C.
*
*   4) Reads in Total Summary Block Table DD if present, as referenced by
*       TotalSummaryOffset in Table C.
*
*   5) Reads in B+ Tree Header Chromosome Index Header Table E, as referenced
*       by chromosomeTreeOffset in Table C.
*
*   6)Reads in B+ Tree Nodes indexing mChromosome ID's for mChromosome regions;
*       Table F for node type (leaf/child), Table G for leaf items,
*       Table H for child node items.
*
*   7) Reads in R+ Tree Chromosome ID Index Header Table K.
*
*   8) Reads in R+ Tree Nodes indexing of data arranged by mChromosome ID's;
*       Table L for node type (leaf or child), Table M for leaf items,
*       Table N for child node items.
*
*   9) Verifies Data Count of data records, as referenced by fullDataOffset in Table C
*
*   10) References data count records of data size defined in Table M of R+ Tree index
*       for all leaf items in the tree.
*
*   11) Reads in zoom level format Table O for each zoom level comprised of
*       zoom record count followed by that many Table P zoom statistics records,
*       followed by an R+ tree of zoom data locations indexed as in Tables L, M, and N.
*
*   12) Returns information on chromosome name keys and chromosome data regions.
*
*   13) Provides iterators using chromosome names and data regions to extract
*       zoom data, Wig data, and Bed data.
*
* */
#include "BBZoomLevelHeader.h"
#include "BBZoomLevelFormat.h"
#include "BBFileHeader.h"
#include "BBZoomLevels.h"
#include "BBTotalSummaryBlock.h"
#include "BPTree.h"
#include "RPTree.h"
#include "BigWigIterator.h"
#include "ZoomLevelIterator.h"



class BBFileReader{

public :

    static const uint64_t BBFILE_HEADER_OFFSET = 0;

    // Defines the Big Binary File (BBFile) access
    private :

    std::string path_;         // BBFile source file/pathname
    std::ifstream* fis_;      // BBFile input stream handle
    int64_t fileOffset_;           // file offset for next item to be read

    BBFileHeader* fileHeader_; // Big Binary file header
    uint32_t dataCount_;             // Number of data records in the file - Table BB
   // bool isLowToHigh;       // BBFile binary data format: low to high or high to low
    uint32_t uncompressBufSize_;     // buffer byte size for data decompression; 0 for uncompressed

    // AutoSQL String defines custom BigBed formats
    uint64_t autoSqlOffset_;
    std::string autoSqlFormat_;

    // This section defines the zoom items if zoom data exists
    uint32_t zoomLevelCount_;       // number of zoom levels defined
    uint64_t zoomLevelOffset_;      // file offset to zoom level headers
    BBZoomLevels* zoomLevels_;   // zoom level headers and data locations

    // Total Summary Block - file statistical info
    uint64_t mTotalSummaryBlockOffset_;
    BBTotalSummaryBlock* totalSummaryBlock_;

    // B+ tree
    uint64_t chromIDTreeOffset_; // file offset to mChromosome index B+ tree
    BPTree* chromosomeIDTree_;     // Container for the mChromosome index B+ tree

    // R+ tree
    uint64_t chromDataTreeOffset_;  // file offset to mChromosome data R+ tree
    RPTree* chromosomeDataTree_;     // Container for the mChromosome data R+ tree



    public:
    BBFileReader();
    virtual ~BBFileReader();

    BBFileReader(std::string path, std::ifstream& stream);
    std::string getBBFilePath();
    std::ifstream* getBBFis();
    BBFileHeader* getBBFileHeader();
    bool isBigBedFile();
     bool isBigWigFile();

    int64_t getChromosomeNameCount();
    int64_t getChromosomeRegionCount();
    int32_t getDecompressionBufSize();


    BBTotalSummaryBlock* getTotalSummaryBlock();
    BPTree* getChromosomeIDTree();
    RPTree* getChromosomeDataTree();
    int32_t getZoomLevelCount();
    BBZoomLevels* getZoomLevels();
    RPChromosomeRegion* getZoomLevelBounds(int32_t zoomLevel, int32_t startChromID, int32_t endChromID);
    RPChromosomeRegion* getZoomLevelBounds(int32_t zoomLevel);
    int32_t getZoomLevelRecordCount(int32_t zoomLevel);
    std::string getChromosomeName(int32_t chromID);
    std::vector<std::string> getChromosomeNames();
    int32_t getChromosomeID(std::string chromosomeKey);
    RPChromosomeRegion* getChromosomeBounds(int32_t startChromID, int32_t endChromID);
    RPChromosomeRegion* getChromosomeBounds();
    std::vector<RPChromosomeRegion*> getChromosomeRegions();
    std::vector<RPChromosomeRegion*> getZoomLevelRegions(int32_t zoomLevel);

    void open(std::string path, std::ifstream& stream);


    ZoomLevelIterator getZoomLevelIterator(unsigned int zoomLevel, std::string startChromosome, int startBase, std::string endChromosome, int endBase, bool contained);
    BigWigIterator getBigWigIterator(std::string startChromosome, int32_t startBase,std::string endChromosome, int32_t endBase, bool contained);
    BigWigIterator getBigWigIterator();
    RPChromosomeRegion* getChromosomeBounds(std::string startChromosome, int startBase,std::string endChromosome, int endBase);
    ZoomLevelIterator getZoomLevelIterator(unsigned int zoomLevel);
    ZoomLevelIterator getZoomLevelIterator(unsigned int zoomLevel, RPChromosomeRegion* selectionRegion, bool contained);

    uint32_t getDataCount();
    uint32_t getDataCount(std::ifstream& fis_, uint64_t fileOffset);

    };


#endif // BBFILEREADER_H_INCLUDED
