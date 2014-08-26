#ifndef BPTREELEAFNODE_H
#define BPTREELEAFNODE_H

// ***************************************************************************
//   BPTreeLeafNode.h (c) 2014
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


#include "BPTreeNode.h"
#include "BPTreeLeafNodeItem.h"
#include <vector>
#include "smart_pointer.h"

class BPTreeLeafNode : public BPTreeNode
{
//    class BPTreeLeafNodeItem;
    public:
        BPTreeLeafNode();
        BPTreeLeafNode(int64_t nodeIndex);
        virtual ~BPTreeLeafNode();

        bool isLeaf() { return isLeafNode_; }
        void setIsLeaf(bool val) { isLeafNode_ = val; }

        bool insertItem(BPTreeNodeItem* item);
        bool deleteItem(int index);

        int32_t getItemCount();
        BPTreeNodeItem* getItem(int index);

        std::vector<BPTreeLeafNodeItem *> getLeafItems();


        int64_t getNodeIndex() { return nodeIndex_; }
        void setNodeIndex(uint64_t val) { nodeIndex_ = val; }


        std::string getLowestChromKey() { return lowestChromKey_; }
        void setLowestChromKey(std::string val) { lowestChromKey_ = val; }
        std::string getHighestChromKey() { return highestChromKey_; }
        void setHighestChromKey(std::string val) { highestChromKey_ = val; }


        int32_t getHighestChromID() { return highestChromID_; }
        void setHighestChromID(int32_t val) { highestChromID_ = val; }
        int32_t getLowestChromID() { return lowestChromID_; }
        void setLowestChromID(int32_t val) { lowestChromID_ = val; }


    protected:
    private:
        bool isLeafNode_;
        uint64_t nodeIndex_;
        std::string lowestChromKey_;
        std::string highestChromKey_;
        uint32_t lowestChromID_;
        uint32_t highestChromID_;
        std::vector<BPTreeLeafNodeItem*> leafItems_;
};

#endif // BPTREELEAFNODE_H
