#ifndef BIGWIGDATABLOCK_H
#define BIGWIGDATABLOCK_H


// ***************************************************************************
//   BigWigDataBlock.h (c) 2014
//   Copyright @ Alexei Nordell-Markovits : Sherbrooke University
//
//    This file is part of the BWReader library.
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

#include "WigItem.h"
#include "RPChromosomeRegion.h"
#include "RPTreeLeafNodeItem.h"
#include <map>
class BigWigDataBlock
{
    public:
        BigWigDataBlock();
        virtual ~BigWigDataBlock();

        std::vector<WigItem> getWigData(RPChromosomeRegion& selectionRegion,std::map<uint32_t, std::string>& chromosomeMap, bool contained);
        BigWigDataBlock(std::ifstream* fis, std::vector<RPTreeLeafNodeItem*>::iterator leafIter, int32_t uncompressBufSize);


    protected:
    private:

    // BigWig data types sizes
    static const int32_t FIXED_STEP_ITEM_SIZE = 4;
    static const int32_t VAR_STEP_ITEM_SIZE = 8;
    static const int32_t BED_GRAPH_ITEM_SIZE = 12;

    // Bed data block access variables   - for reading in bed records from a file
    int64_t fileOffset_;       // Wig data block file offset
    int64_t leafDataSize_;     // byte size for data block specified in the R+ leaf

    // defines the bigWig data source
  //  std::map<uint32_t, std::string>&  chromosomeMap_;  // map of chromosome ID's and corresponding names
    std::vector<RPTreeLeafNodeItem*>::iterator leafHitItem_;   // R+ leaf item containing data block location

    // uncompressed byte stream buffer and readers
  //  char* wigBuffer_;      // buffer containing leaf block data uncompressed
    std::vector<char> wigBuffer_;
    int32_t remDataSize_;       // number of uncompressed data bytes not extracted

    // Wig data extraction members
    std::vector<WigItem> wigItemList_;  // array of Wig section items

};

#endif // BIGWIGDATABLOCK_H
