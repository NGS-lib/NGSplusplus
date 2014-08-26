
// ***************************************************************************
//   BigWigSection.cpp (c) 2014
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

#include "BigWigSection.h"
#include "BigWigSectionHeader.h"
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include "endian_helper.h"
#include <stdlib.h>
BigWigSection::BigWigSection()
{
    //ctor
}

BigWigSection::~BigWigSection()
{

    //TODO free memory here
//for(int i = 0; i < wigItemList_.size(); ++i)
//    delete wigItemList[i];

}


 /*
    *   Constructor for a BigWig data section which includes the section header
    *   and Wig data items.
    *
    *   Parameters:
    *       sectionBuffer - buffer contains decompressed Wig section header + data
    *       sectionIndex - wig section index for leaf data block
    *       chromIDTree - B+ chromosome index tree returns chromosome names for ID's
    *       isLowToHigh - if true, data byte order is low to high ; else is high to low
    *       leafHitItem - contains leaf node information for testing against selection region
    *
    * */
    BigWigSection::BigWigSection(std::vector<char>& sectionBuffer, std::map<uint32_t,std::string>* chromosomeMap,
                         std::vector<RPTreeLeafNodeItem*>::iterator leafHitItem):leafHitItem_(leafHitItem),chromosomeMap_(chromosomeMap),sectionBuffer_(sectionBuffer){

//    typedef WigTypeNamespace::WigItemType WigItemType;


    try {

        myStream_.rdbuf()->pubsetbuf(&sectionBuffer[0],sectionBuffer.size());

        // wrap the Wig section buffer as an input stream and get the section header
        // Note: A RuntimeException is thrown if header is not read properly

        wigSectionHeader_ = new BigWigSectionHeader(myStream_);


        // check for valid Wig item type
        if(wigSectionHeader_->getItemType() == WigTypeNamespace::Unknown)
        {
            throw std::runtime_error("Read error on wig section leaf index ");
        }

        // include header in data segment size accounting
        sectionDataSize_ = wigSectionHeader_->SECTION_HEADER_SIZE;

        // use method getSectionData to extract section data

        }
        catch(...)
        {
            throw;
        }
    }

    /*
    *   Method returns the if the Wig items defined in this section are valid.
    *
    *   Note: Use BigWigSectionHeader to obtain more information on
    *       Wig section data specifications.
    *
    *   Returns:
    *       Specifies if Wig section has a valid data item type.
    * */
    bool BigWigSection::isValidSectionType(){
        return wigSectionHeader_->IsValidType();
    }

    /*
    *   Method returns the Wig Section Header
    *
    *   Returns:
    *       Wig section header
    * */
    int32_t BigWigSection::getItemCount() {
        return wigSectionHeader_->getItemCount();
    }

    /*
    *   Method returns the Wig Section Header
    *
    *   Returns:
    *       Wig section header
    * */
    BigWigSectionHeader* BigWigSection::getSectionHeader() {
        return wigSectionHeader_;
    }

    /*
    *   Method returns the number bytes of decompressed data in this section.
    *
    *   Returns:
    *       Number of uncompressed bytes read for the Wig data section
    * */
    int32_t BigWigSection::BigWigSection::getSectionDataSize() {
        return sectionDataSize_;
    }

    /*
    *   Method reads Wig data items within the decompressed block buffer for the selection region.
    *
    *   Parameters:
    *       selectionRegion - chromosome selection region for item extraction
    *       contained - indicates select region must be contained in value region
    *           if true, else may intersect selection region for extraction
    *
    *   Returns:
    *     Size in bytes for the wig data section.
    *     Items read in the wig segment data block are added to the wig item list .
    *
     *   Note: Unlike ZoomLevel and BigBed formats, the Wig Section data block header contains
    *   an item count used to determine the end of data read.
    * */
    int32_t BigWigSection::getSectionData(RPChromosomeRegion* selectionRegion, bool contained,
                              std::vector<WigItem*>& wigItemList) {

        // get the section's data item specifications
        // Note: A RuntimeException is thrown if wig section header is not read properly
        int32_t chromID =  wigSectionHeader_->getChromID();
        std::string chromosome = (*chromosomeMap_)[chromID];
        int32_t itemCount = wigSectionHeader_->getItemCount();
        int32_t chromStart = wigSectionHeader_->getChromosomeStart();
        int32_t chromEnd = wigSectionHeader_->getChromosomeEnd();
        int32_t itemStep = wigSectionHeader_->getItemStep();
        int32_t itemSpan =  wigSectionHeader_->getItemSpan();
        int32_t itemIndex = 0;
        int32_t startBase = 0;
        int32_t endBase = 0;
        float value = 0.0f;

        // find Wig data type - BBFile Table J item type
        WigTypeNamespace::WigItemType itemType = wigSectionHeader_->getItemType();

        // check if all leaf items are selection hits
        RPChromosomeRegion* itemRegion = new RPChromosomeRegion(chromID, chromStart,
                            chromID, chromEnd);
        int32_t leafHitValue = itemRegion->compareRegions(selectionRegion);


        // extract Wig data records
        // Note: the buffer input stream is positioned past section header
        try {
            for(int32_t index = 0; index < itemCount; ++index) {
                ++itemIndex;


                    if(itemType == WigTypeNamespace::FixedStep){
                        startBase = chromStart;
                        endBase = startBase + itemSpan;


                        myStream_.read( reinterpret_cast<char*>(&value) , sizeof(float) );
                        value=endian::LittleFloat(value);
                        chromStart = startBase + itemStep;
                        sectionDataSize_ += BigWigSectionHeader::FIXEDSTEP_ITEM_SIZE;
                    }
                    else if(itemType == WigTypeNamespace::VarStep){


                        myStream_.read( reinterpret_cast<char*>(&startBase) , sizeof(uint32_t) );
                        startBase=endian::LittleLong(startBase);

                        endBase = startBase + itemSpan;

                        myStream_.read( reinterpret_cast<char*>(&value) , sizeof(float) );
                        value=endian::LittleFloat(value);
                        sectionDataSize_ += BigWigSectionHeader::VARSTEP_ITEM_SIZE;


                    }
                    else if(itemType == WigTypeNamespace::BedGraph){

                        myStream_.read( reinterpret_cast<char*>(&startBase) , sizeof(uint32_t) );
                        startBase=endian::LittleLong(startBase);

                        myStream_.read( reinterpret_cast<char*>(&endBase) , sizeof(uint32_t) );
                        endBase=endian::LittleLong(endBase);

                        myStream_.read( reinterpret_cast<char*>(&value) , sizeof(float) );
                        value=endian::LittleFloat(value);

                        sectionDataSize_ += BigWigSectionHeader::BEDGRAPH_ITEM_SIZE;
                    }

                // contained leaf region items are always added - otherwise test conditions
                if(leafHitValue == 0) {
                    WigItem* bbItem = new WigItem(itemIndex, chromosome, startBase, endBase, value);
                    wigItemList.push_back(bbItem);
                }
                else {
                    itemRegion = new RPChromosomeRegion(chromID, startBase, chromID, endBase);
                    int32_t itemHitValue = itemRegion->compareRegions(selectionRegion);

                    // hitValue < 2 needed for intersection; hitValue < 1 needed for contained = true
                    if(itemHitValue == 0 || ( !contained && abs(itemHitValue) < 2) ) {
                        WigItem* bbItem = new WigItem(itemIndex, chromosome, startBase, endBase, value);
                        wigItemList.push_back(bbItem);
                    }
                }

            }

        }catch(...) {
            throw new std::runtime_error("Read error for Wig section item ");
        }

        return sectionDataSize_;
    }


