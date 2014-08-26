#ifndef RPTreeHeader_H
#define RPTreeHeader_H


// ***************************************************************************
//   RPTreeHeader.h (c) 2014
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



#include "stdint.h"
#include <iostream>
#include <fstream>
#include "endian_helper.h"


class RPTreeHeader
{
    public:

        static const uint32_t RPTREE_HEADER_SIZE = 48;
        static const uint32_t RPTREE_MAGIC_LTH = 0x2468ACE0;
        static const uint32_t RPTREE_MAGIC_HTL = 0xE0AC6824;



        RPTreeHeader();
        virtual ~RPTreeHeader();
        RPTreeHeader(std::ifstream& fis, uint64_t fileOffset);
        bool readHeader(std::ifstream& fis, uint64_t fileOffset);


        uint64_t getRPTreeOffset_() { return rpTreeOffset_; }
        void setRPTreeOffset(uint64_t val) { rpTreeOffset_ = val; }
        bool isHeaderOK() { return headerOK_; }


        uint32_t getMagic() { return magic_; }
        void setMagic(uint32_t val) { magic_ = val; }
        uint32_t getBlockSize() { return blockSize_; }
        void setBlockSize(uint32_t val) { blockSize_ = val; }
        uint64_t getItemCount() { return itemcount_; }


        uint32_t getStartChromID() { return startChromID_; }
        void setStartChromID(uint32_t val) { startChromID_ = val; }
        uint32_t getStartBase() { return startBase_; }
        void setStartBase(uint32_t val) { startBase_ = val; }
        uint32_t getEndChromID() { return endChromID_; }
        void setEndChromID(uint32_t val) { endChromID_ = val; }
        uint32_t getEndBase() { return endBase_; }
        void setEndBase(uint32_t val) { endBase_ = val; }
        uint32_t getEndFileOffset() { return endFileOffset_; }
        void setEndFileOffset(uint32_t val) { endFileOffset_ = val; }
        uint32_t getItemsPerSlot() { return itemsPerSlot_; }
        void setItemsPerSlot(uint32_t val) { itemsPerSlot_ = val; }
        uint32_t getReserved() { return reserved_; }
        void setReserved(uint32_t val) { reserved_ = val; }


        uint32_t getHeaderSize() {
            return RPTREE_HEADER_SIZE;
        }
    protected:
    private:


        uint64_t rpTreeOffset_;
        bool headerOK_;
        uint32_t magic_;
        uint32_t blockSize_;
        uint64_t itemcount_;
        uint32_t startChromID_;
        uint32_t startBase_;
        uint32_t endChromID_;
        uint32_t endBase_;
        uint64_t endFileOffset_;
        uint32_t itemsPerSlot_;
        uint32_t reserved_;
};

#endif // RPTreeHeader_H
