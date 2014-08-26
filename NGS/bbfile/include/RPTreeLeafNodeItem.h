#ifndef RPTREELEAFNODEITEM_H
#define RPTREELEAFNODEITEM_H


// ***************************************************************************
//   RPTreeLeafNodeItem.h (c) 2014
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
#include "RPTreeNodeItem.h"
#include "RPChromosomeRegion.h"

class RPTreeLeafNodeItem : public RPTreeNodeItem, public RPChromosomeRegion
{
    public:
        //RpTreeLeafNodeItem();

        using RPChromosomeRegion::compareRegions;



        RPTreeLeafNodeItem(uint32_t startChromID, uint32_t startBase,
                              uint32_t endChromID, uint32_t endBase, uint64_t dataOffset, uint64_t dataSize);
        virtual ~RPTreeLeafNodeItem();
        uint64_t getDataOffset() { return dataOffset_; }
    //    void setDataOffset(uint64_t val) { dataOffset_ = val; }
        uint64_t getDataSize() { return dataSize_; }
    //    void setDataSize(uint64_t val) { dataSize_ = val; }


        RPChromosomeRegion* getChromosomeBounds();

    protected:
    private:
        uint64_t dataOffset_;
        uint64_t dataSize_;
};

#endif // RPTREELEAFNODEITEM_H
