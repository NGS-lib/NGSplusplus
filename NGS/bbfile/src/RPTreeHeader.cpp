// ***************************************************************************
//   RpTreeHeader.cpp (c) 2014
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

#include "RPTreeHeader.h"
#include "endian_helper.h"
#include <stdexcept>

RPTreeHeader::RPTreeHeader()
{
    //ctor
}

RPTreeHeader::~RPTreeHeader()
{
    //dtor
}

 // constructor   - reads from file input stream
    /*
    *   Constructor
    *
    *   Parameters:
    *       fis - file input stream handle
    *       fileOffset - file offset to the RP tree header
    *       isLowToHigh - if true, indicates low to high byte order, else high to low
    * */
RPTreeHeader::RPTreeHeader(std::ifstream& fis, uint64_t fileOffset) {

//    uint64_t itemsCount;

   rpTreeOffset_ =  fileOffset;

   // Note: a bad R+ Tree header will result in false returned
   headerOK_ =  readHeader(fis, rpTreeOffset_);
}

  /*
  * Reads in the R+ Tree Header.
  *
  * Returns status of for tree header read; true if read, false if not.
  * */
   bool RPTreeHeader::readHeader(std::ifstream& fis, uint64_t fileOffset){

       try {
           // Read R+ tree header into a buffer

        fis.clear();
        fis.seekg(fileOffset);
        if (fis.eof())
                std::cerr<<"Hit end of file in seekg in RPTreeHeader::readHeader\n";


        fis.read( reinterpret_cast<char*>(&magic_) , sizeof(uint32_t) );
        magic_=endian::LittleLong(magic_);

       // check for a valid B+ Tree Header
       if(magic_ != RPTREE_MAGIC_LTH)
           return false;

       // Get mChromosome B+ header information
        fis.read( reinterpret_cast<char*>(&blockSize_) , sizeof(uint32_t) );
        blockSize_=endian::LittleLong(blockSize_);

        fis.read( reinterpret_cast<char*>(&itemcount_) , sizeof(uint64_t) );
        itemcount_=endian::LittleDouble(itemcount_);

        fis.read( reinterpret_cast<char*>(&startChromID_) , sizeof(uint32_t) );
        startChromID_=endian::LittleLong(startChromID_);

        fis.read( reinterpret_cast<char*>(&startBase_) , sizeof(uint32_t) );
        startBase_=endian::LittleLong(startBase_);

        fis.read( reinterpret_cast<char*>(&endChromID_) , sizeof(uint32_t) );
        endChromID_=endian::LittleLong(endChromID_);

        fis.read( reinterpret_cast<char*>(&endBase_) , sizeof(uint32_t) );
        endBase_=endian::LittleLong(endBase_);


        fis.read( reinterpret_cast<char*>(&endFileOffset_) , sizeof(uint64_t) );
        endFileOffset_=endian::LittleDouble(endFileOffset_);

        fis.read( reinterpret_cast<char*>(&itemsPerSlot_) , sizeof(uint32_t) );
        itemsPerSlot_=endian::LittleLong(itemsPerSlot_);

        fis.read( reinterpret_cast<char*>(&reserved_) , sizeof(uint32_t) );
        reserved_=endian::LittleLong(reserved_);



       }
       catch(...) {
               throw new std::runtime_error("Error reading R+ tree header ");
           }

       // success
        return true;
   }
