// ***************************************************************************
//   RPTreeLeafNodeItem.cpp (c) 2014
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

#include "RPTreeLeafNodeItem.h"

RPTreeLeafNodeItem::~RPTreeLeafNodeItem()
{
    //dtor
}


 /*  Constructor for leaf node items.
    *
    *   Parameters:
    *       itemIndex - index of item belonging to a leaf node
    *       startChromID - starting chromosome/contig for item
    *       startBase - starting base for item
    *       endChromID - ending chromosome/contig for item
    *       endBase - ending base for item
    *       dataOffset - file location for leaf chromosome/contig data
    *       dataSize - size of (compressed) leaf data region in bytes
    *
    * */
    //TODO This was not assigning chrom values, was this the problem?
    RPTreeLeafNodeItem::RPTreeLeafNodeItem(uint32_t startChromID, uint32_t startBase, uint32_t endChromID, uint32_t endBase, uint64_t dataOffset, uint64_t dataSize):
        RPChromosomeRegion(startChromID,startBase,endChromID,endBase), dataOffset_(dataOffset),dataSize_(dataSize)
    {



    }

    RPChromosomeRegion* RPTreeLeafNodeItem::getChromosomeBounds() {
        return this;
    }


