
// ***************************************************************************
//   BigWigSectionHeader.cpp (c) 2014
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

#include "BigWigSectionHeader.h"
#include <fstream>
#include <iostream>
#include "endian_helper.h"
#include <stdexcept>
#include <sstream>

BigWigSectionHeader::BigWigSectionHeader()
{
    //ctor
}

BigWigSectionHeader::~BigWigSectionHeader()
{
    //dtor
}

    /*
    *   Constructor creates a Wig Section Header (Table J) from uncompressed buffer.
    *
    *   Parameters:
    *       mLbdis - // buffer stream containing section header arranged high to low bytes
    * */
    BigWigSectionHeader::BigWigSectionHeader(std::ifstream& bdis) {
        // get Wig Section Header
       char type;
        try {

            bdis.read( reinterpret_cast<char*>(&chromID_) , sizeof(uint32_t) );
            chromID_=endian::LittleLong(chromID_);

            bdis.read( reinterpret_cast<char*>(&chromStart_) , sizeof(uint32_t) );
            chromStart_=endian::LittleLong(chromStart_);

            bdis.read( reinterpret_cast<char*>(&chromEnd_) , sizeof(uint32_t) );
            chromEnd_=endian::LittleLong(chromEnd_);

            bdis.read( reinterpret_cast<char*>(&itemStep_) , sizeof(uint32_t) );
            itemStep_=endian::LittleLong(itemStep_);

            bdis.read( reinterpret_cast<char*>(&itemSpan_) , sizeof(uint32_t) );
            itemSpan_=endian::LittleLong(itemSpan_);

            bdis.read( reinterpret_cast<char*>(&type) , sizeof(unsigned char) );
            type=endian::LittleByte(type);

            bdis.read( reinterpret_cast<char*>(&reserved_) , sizeof(unsigned char) );
            reserved_=endian::LittleByte(reserved_);

            bdis.read( reinterpret_cast<char*>(&itemCount_) , sizeof(uint16_t) );
            itemCount_=endian::LittleShort(itemCount_);

        }catch(...) {
            throw std::runtime_error("Error reading wig section header");
        }

        // tag as valid
        isValidType_ = getItemType(type);
    }

    BigWigSectionHeader::BigWigSectionHeader(std::stringstream& bdis){

 // get Wig Section Header

       unsigned char type;
        try {

            bdis.read( reinterpret_cast<char*>(&chromID_) , sizeof(uint32_t) );
            chromID_=endian::LittleLong(chromID_);


            bdis.read( reinterpret_cast<char*>(&chromStart_) , sizeof(uint32_t) );
            chromStart_=endian::LittleLong(chromStart_);

            bdis.read( reinterpret_cast<char*>(&chromEnd_) , sizeof(uint32_t) );
            chromEnd_=endian::LittleLong(chromEnd_);

            bdis.read( reinterpret_cast<char*>(&itemStep_) , sizeof(uint32_t) );
            itemStep_=endian::LittleLong(itemStep_);

            bdis.read( reinterpret_cast<char*>(&itemSpan_) , sizeof(uint32_t) );
            itemSpan_=endian::LittleLong(itemSpan_);

            bdis.read( reinterpret_cast<char*>(&type) , sizeof(unsigned char) );
            type=endian::LittleByte(type);

            bdis.read( reinterpret_cast<char*>(&reserved_) , sizeof(unsigned char) );
            reserved_=endian::LittleByte(reserved_);

            bdis.read( reinterpret_cast<char*>(&itemCount_) , sizeof(uint16_t) );
            itemCount_=endian::LittleShort(itemCount_);

        }catch(...) {
            throw std::runtime_error("Error reading wig section header");
        }

        // tag as valid
        isValidType_ = getItemType(type);



    }




    /*
    *   Method returns the chromosome ID
    *
    *   Returns:
    *       Chromosome ID for the section's region
    * */
    uint32_t BigWigSectionHeader::getChromID() {
        return chromID_;
    }

    /*
    *   Method returns the chromosome starting base
    *
    *   Returns:
    *       Chromosome start base for the section's region
    * */
    uint32_t BigWigSectionHeader::getChromosomeStart() {
        return chromStart_;
    }

    /*
    *   Method returns the chromosome ending base
    *
    *   Returns:
    *       Chromosome end base for the section's region
    * */
    uint32_t BigWigSectionHeader::getChromosomeEnd() {
        return chromEnd_;
    }

    /*
    *   Method returns the base pairs step between items.
    *
    *   Returns:
    *       Chromosome base step between fixed step sections
    * */
    uint32_t BigWigSectionHeader::getItemStep() {
        return itemStep_;
    }

    /*
    *   Method returns the base pairs span in items.
    *
    *   Returns:
    *       Chromosome base span for fixed and variable step sections
    * */
    uint32_t BigWigSectionHeader::getItemSpan() {
        return itemSpan_;
    }

    /*
    *   Method returns the item type for the section's Wig data.
    *
    *   Returns:
    *       Section item type for Wig data
    * */
    WigTypeNamespace::WigItemType BigWigSectionHeader::getItemType() {
        return itemType_;
    }

    /*
    *   Method returns if the section's data item type is valid.
    *
    *   Returns:
    *       Specifies if section's data iytem type is valid
    * */
    bool BigWigSectionHeader::IsValidType() {
        return isValidType_;
    }

    /*
    *   Method returns the number of section items.
    *
    *   Returns:
    *       Number of items defined for the section
    * */
    uint16_t BigWigSectionHeader::getItemCount() {
        return itemCount_;
    }

    /*
    *   Method returns the reserved value for the section.
    *
    *   Returns:
    *       Reserved byte for the section (should always be 0)
    * */
    unsigned char BigWigSectionHeader::getReserved() {
        return reserved_;
    }

    /*
    *   Method determines the Wig data type.
    *
    *   Parameters:
    *       byte type read from Wig section header
    *
    *   Returns:
    *       Indicates if type is a valid Wig item type
    * */
    bool BigWigSectionHeader::getItemType(char type){
        bool isValid;

        if(type == 1){
            itemType_ = WigTypeNamespace::BedGraph;
            itemDescription_ = "Wig Bed Graph";
            isValid = true;
        }
        else if(type == 2){
            itemType_ = WigTypeNamespace::VarStep;
            itemDescription_ = "Wig Variable Step";
            isValid = true;
        }
        else if(type == 3){
            itemType_ = WigTypeNamespace::FixedStep;
            itemDescription_ = "Wig Fixed Step";
            isValid = true;
        }
        else {
            itemType_ = WigTypeNamespace::Unknown;
            itemDescription_ = "Wig Type Unknown";
            isValid = false;
        }

        return isValid;
    }

