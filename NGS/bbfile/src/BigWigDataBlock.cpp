
// ***************************************************************************
//   BigWigDataBlock.cpp (c) 2014
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

#include "BigWigDataBlock.h"
#include "BigWigSection.h"
#include "BigWigSectionHelper.h"
#include "decompress_util.h"
#include <vector>
#include <string>


BigWigDataBlock::BigWigDataBlock()
{
    //ctor
}

BigWigDataBlock::~BigWigDataBlock()
{
 // for(uint32_t i = 0; i < wigItemList_.size(); ++i)
 //   delete wigItemList_[i];
}

//std::vector<RPTreeLeafNodeItem>::iterator
  /*
    *   Constructor for Wig data block reader.
    *
    *   Parameters:
    *       fis - file input stream handle
    *       leafHitItem - R+ tree leaf hit item containing data block file location and hit status
    *       chromIDTree - B+ chromosome index tree returns chromosome ID's for names
    *       uncompressBufSize - byte size for decompression buffer; else 0 for uncompressed
    *
    * */
    BigWigDataBlock::BigWigDataBlock(std::ifstream* fis, std::vector<RPTreeLeafNodeItem*>::iterator leafIter, int32_t uncompressBufSize){

    //TODO this is in reafctoring/test
        this->leafHitItem_ = leafIter;
      //  this->chromosomeMap_ = &chromosomeMap;

        RPTreeLeafNodeItem* itrTest =*leafHitItem_;

        fileOffset_ = itrTest->getDataOffset();
        leafDataSize_ = itrTest->getDataSize();

        fis->clear();
        fis->seekg(fileOffset_);

        try {
            wigBuffer_.clear();
            wigBuffer_.resize(uncompressBufSize);

            return_compressed(*fis,&wigBuffer_[0],uncompressBufSize);
           //read_compressed(*fis,&wigBuffer_[0],leafDataSize_);

/***********************************************
            // decompress if necessary - the buffer size is 0 for uncompressed data
            // Note:  BBFile Table C specifies a decompression buffer size
********************************************/
        }catch(...) {
             throw std::runtime_error("Error reading Wig section for leaf item");
        }
        // initialize unread data size
        remDataSize_ = leafDataSize_;
    }

    /*
    *   Method reads all Wig data sections within the decompressed block buffer
    *   and returns those items in the chromosome selection region.
    *
    *   Parameters:
    *       selectionRegion - chromosome region for selecting Wig values
   *       contained - indicates selected data must be contained in selection region
    *           if true, else may intersect selection region
    *
    *   Returns:
    *      Wig sections in selected from the data block; else null for none selected.
    *
    * */
    std::vector<WigItem> BigWigDataBlock::getWigData(RPChromosomeRegion& selectionRegion,std::map<uint32_t, std::string>& chromosomeMap,
                                                       bool contained)
  {
    try {
        for(int index = 0; remDataSize_ > 0; ++index) {
            // extract items in the Wig data section
            // Note: A RuntimeException is thrown if wig section is not read properly

            //TODO: this is Much to Java-like. There is no reason to create a "throaway object" to decompress the data
            //when we could call a free function.

            int sectionBytes= extractSectionData(selectionRegion,wigBuffer_,chromosomeMap,contained,wigItemList_);


         //   BigWigSection* wigSection = new BigWigSection(wigBuffer_, chromosomeMap_, leafHitItem_);

            // get wig section items and section bytes read
          //  int sectionBytes = wigSection->getSectionData(selectionRegion, contained, wigItemList_);

            //TODO remove this OMGWTF delete once we refactor
         //   delete (wigSection);
            // adjust remaining data block size

            remDataSize_ -= sectionBytes;
        }

        return wigItemList_;
    }catch(...){
    throw;

    }

    }


