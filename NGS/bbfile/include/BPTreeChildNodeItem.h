#ifndef BPTREECHILDNODEITEM_H
#define BPTREECHILDNODEITEM_H

// ***************************************************************************
//   BPTreeChildNodeItem.h (c) 2014
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
#include "BPTreeNode.h"

class BPTreeChildNodeItem : public BPTreeNodeItem
{
    public:
     //   BPTreeChildNodeItem();
        virtual ~BPTreeChildNodeItem();
        BPTreeChildNodeItem(uint32_t itemIndex, std::string chromKey, BPTreeNode* childNode);

        std::string getChromKey() { return chromKey_; }
        void setChromKey(std::string val) { chromKey_ = val; }
        BPTreeNode* getChildNode() { return childNode_; }
        void setChildNode(BPTreeNode* val) { childNode_ = val; }
        uint64_t getItemIndex() { return itemIndex_; }
        void setItemIndex(uint64_t val) { itemIndex_ = val; }
        bool isLeafItem() { return isLeafItem_; }
        void setIsLeafItem(bool val) { isLeafItem_ = val; }
        bool chromKeysMatch(std::string chromKey) ;
    protected:
    private:
        std::string chromKey_;
        BPTreeNode* childNode_;
        uint64_t itemIndex_;
        bool isLeafItem_;
};

#endif // BPTREECHILDNODEITEM_H
