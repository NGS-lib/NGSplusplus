#ifndef BPTREELEAFNODEITEM_H
#define BPTREELEAFNODEITEM_H


// ***************************************************************************
//   BPTreeLeafNodeItem.h (c) 2014
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


#include "BPTreeNodeItem.h"
#include "BPTreeLeafNodeItem.h"

class BPTreeLeafNodeItem : public BPTreeNodeItem
{
    public:
        BPTreeLeafNodeItem();
        virtual ~BPTreeLeafNodeItem();

        BPTreeLeafNodeItem(uint64_t leafIndex, std::string chromKey, uint32_t chromID, uint32_t chromSize);


        bool chromKeysMatch(std::string chromKey);

        bool isLeafItem() { return isLeafItem_; }
        void setIsLeafItem(bool val) { isLeafItem_ = val; }
        uint64_t getLeafIndex() { return leafIndex_; }
        void setLeafIndex(uint64_t val) { leafIndex_ = val; }
        std::string getChromKey() { return chromKey_; }
        void setChromKey(std::string val) { chromKey_ = val; }

        uint32_t getChromID() { return chromID_; }
        void setChromID(uint32_t val) { chromID_ = val; }


        uint32_t setChromSize() { return chromSize_; }
        void setChromSize(uint32_t val) { chromSize_ = val; }

        uint64_t getItemIndex(){return leafIndex_;}

    protected:
    private:
        bool isLeafItem_;
        uint64_t leafIndex_;
        std::string chromKey_;
        uint32_t chromID_;
        uint32_t chromSize_;


};

#endif // BPTREELEAFNODEITEM_H