#ifndef BPTREEHEADER_H
#define BPTREEHEADER_H

// ***************************************************************************
//   BPTreeHeader.h (c) 2014
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


#include <string>
#include <iostream>
#include "stdint.h"
class BPTreeHeader
{
    public:

        static const uint32_t BPTREE_HEADER_SIZE = 32;
        static const uint32_t BPTREE_MAGIC_LTH = 0x78CA8C91;

        BPTreeHeader(std::ifstream& fis, uint64_t fileOffset);
        uint64_t GetheaderOffset_() { return headerOffset_; }
        void setHeaderOffset(uint64_t val) { headerOffset_ = val; }
        bool isHeaderOK() { return headerOk_; }
        void SetHeaderOk(bool val) { headerOk_ = val; }
        uint32_t getMagic() { return magic_; }
        void setMagic(uint32_t val) { magic_ = val; }
        uint32_t getKeySize() { return keySize_; }
        void setKeySize(uint32_t val) { keySize_ = val; }
        uint32_t getBlockSize() { return blockSize_; }
        void setBlockSize(uint32_t val) { blockSize_ = val; }
        uint32_t getValSize() { return valSize_; }
        void setValSize(uint32_t val) { valSize_ = val; }
        uint64_t getItemCount() { return itemCount_; }
        void setItemCount(uint64_t val) { itemCount_ = val; }
        uint64_t getReserved() { return reserved_; }
        void setReserved(uint64_t val) { reserved_ = val; }
    protected:
    private:

        bool readHeader(std::ifstream& fis, uint64_t fileOffset);

        uint64_t headerOffset_;
        bool headerOk_;
        uint32_t magic_;
        uint32_t keySize_;
        uint32_t blockSize_;
        uint32_t valSize_;
        uint64_t itemCount_;
        uint64_t reserved_;
};

#endif // BPTREEHEADER_H
