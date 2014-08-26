// ***************************************************************************
//   BPTreeChildNodeItem.cpp (c) 2014
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

#include "BPTreeChildNodeItem.h"

//BPTreeChildNodeItem::BPTreeChildNodeItem()
//{
//    isLeafItem_=false;
//}

BPTreeChildNodeItem::~BPTreeChildNodeItem()
{
    //dtor
}


 /*
    *   Constructs a B+ tree child node item with the supplied information.
    *
    *   Parameters:
    *       itemIndex - node item index
    *       chromKey - chromosome name key
    *       childNode - assigned child node object
    * */
    BPTreeChildNodeItem::BPTreeChildNodeItem(uint32_t itemIndex, std::string chromKey, BPTreeNode* childNode){
        isLeafItem_= false;
        itemIndex_ = itemIndex;
        chromKey_  = chromKey;
        childNode_ = childNode;
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
    bool BPTreeChildNodeItem::chromKeysMatch(std::string chromKey) {
        std::string thisKey = chromKey_;
        std::string thatKey = chromKey;

        // Note: must have the same length to compare chromosome names
        uint32_t thisKeyLength = thisKey.size();
        uint32_t thatKeyLength = thatKey.size();

        // check if need to truncate the larger string
        if(thisKeyLength > thatKeyLength)
            thisKey = thisKey.substr(0,thatKeyLength);
        else if(thatKeyLength > thisKeyLength)
            thatKey = thatKey.substr(0,thisKeyLength);

        if (thisKey==thatKey)
            return true;
        else
            return false;
    }


