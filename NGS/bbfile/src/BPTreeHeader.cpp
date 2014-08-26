// ***************************************************************************
//   BPTreeHeader.cpp (c) 2014
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

#include "BPTreeHeader.h"
#include <iostream>
#include <fstream>
#include "stdint.h"
#include <stdexcept>
#include "endian_helper.h"
 /*
   *    Constructor for reading in a B+ tree header a from a file input stream.
   *
   *    Parameters:
   *        fis - file input handle
   *        fileOffset - file offset to the B+ tree header
   *        isLowToHigh - indicates byte order is low to high, else is high to low
   * */
    BPTreeHeader::BPTreeHeader(std::ifstream& fis, uint64_t fileOffset) {

       // save the seekable file handle  and B+ Tree file offset
       headerOffset_ = fileOffset;

       // Note: a bad B+ Tree header will result in false returned
       headerOk_ =  readHeader(fis, headerOffset_);
    }

   /*
   * Reads in the B+ Tree Header.
   * Returns status of B+ tree header read; true if read, false if not.
   * */
    bool BPTreeHeader::readHeader(std::ifstream& fis, uint64_t fileOffset) {

        try {
            // Read B+ tree header into a buffer
            fis.clear();
            fis.seekg(fileOffset, std::ios_base::beg);
            if (fis.eof())
                std::cerr<<"Hit end of file in seekg in BPTreeHeader::readHeader\n";
            // decode header
            fis.read( reinterpret_cast<char*>(&magic_) , sizeof(uint32_t) );
            magic_=endian::LittleLong(magic_);

            if(magic_ != BPTREE_MAGIC_LTH)
                return false;
                // Get mChromosome B+ header information


            fis.read( reinterpret_cast<char*>(&blockSize_) , sizeof(uint32_t) );
            blockSize_=endian::LittleLong(blockSize_);

            fis.read( reinterpret_cast<char*>(&keySize_) , sizeof(uint32_t) );
            keySize_=endian::LittleLong(keySize_);


            fis.read( reinterpret_cast<char*>(&valSize_) , sizeof(uint32_t) );
            valSize_=endian::LittleLong(valSize_);


            fis.read( reinterpret_cast<char*>(&itemCount_) , sizeof(uint64_t) );
            itemCount_=endian::LittleDouble(itemCount_);


            fis.read( reinterpret_cast<char*>(&reserved_) , sizeof(uint64_t) );
            reserved_=endian::LittleDouble(reserved_);


        }catch(...) {
            throw std::runtime_error("Error reading B+ tree header \n");
            }

         return true;
    }
