// ***************************************************************************
//   BPTreeLeafNodeItem.cpp (c) 2014
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

#include "BPTreeLeafNodeItem.h"
#include <string>
#include "stdint.h"
#include <algorithm>
BPTreeLeafNodeItem::BPTreeLeafNodeItem()
{
    //ctor
    isLeafItem_=true;
}

BPTreeLeafNodeItem::~BPTreeLeafNodeItem()
{
    //dtor
}


  /*
    *   Constructs a B+ tree leaf node item with the supplied information.
    *
    *   Parameters:
    *       leafIndex - leaf item index
    *       chromKey - chromosome/contig name key
    *       chromID - chromosome ID assigned to the chromosome name key
    *       chromsize - number of bases in the chromosome/contig
    * */
     BPTreeLeafNodeItem::BPTreeLeafNodeItem(uint64_t leafIndex, std::string chromKey, uint32_t chromID, uint32_t chromSize) {

         isLeafItem_=true;
        leafIndex_ = leafIndex;
        chromKey_ = chromKey;
        chromID_ = chromID;
        chromSize_ = chromSize;
    }


    /*
   *   Method compares supplied chromosome key with leaf node key.
   *
   *   Parameters:
   *       chromKey - chromosome name ley to compare
   *
   *   Returns:
   *       true, if keys are equal; false if keys are different
   * */
    bool BPTreeLeafNodeItem::chromKeysMatch(std::string chromKey) {
        std::string thisKey = this->chromKey_;
        //Remove whitespaces

        std::string::iterator end_pos = std::remove(thisKey.begin(), thisKey.end(), '\000');
        thisKey.erase(end_pos, thisKey.end());

        return (thisKey==chromKey);
    }


