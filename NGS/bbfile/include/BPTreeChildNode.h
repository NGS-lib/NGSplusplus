#ifndef BPTREECHILDNODE_H
#define BPTREECHILDNODE_H

// ***************************************************************************
//   BPTreeChildNode.h (c) 2014
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
#include <vector>
#include "BPTreeChildNodeItem.h"
class BPTreeChildNode : public BPTreeNode
{
    public:
        BPTreeChildNode();
        virtual ~BPTreeChildNode();

        BPTreeChildNode(int64_t nodeIndex);

        int64_t getNodeIndex() { return nodeIndex_; };
        void setNodeIndex(int64_t val) { nodeIndex_ = val; };

        bool isLeaf(){return isLeafNode_;};

        bool insertItem(BPTreeNodeItem* item);
        bool deleteItem(int index);
        BPTreeNodeItem* getItem(int index);

        int32_t getItemCount();

        std::string getLowestChromKey();
        std::string getHighestChromKey();
        int32_t getLowestChromID();
        int32_t getHighestChromID();

          std::vector<BPTreeChildNodeItem*> getChildItems();

    protected:
    private:

        bool isLeafNode_;
        int64_t nodeIndex_;
        int32_t lowestChromID_;
        int32_t highestChromID_;
        std::string highestChromKey_;
        std::string lowestChromKey_;
        std::vector<BPTreeChildNodeItem *> childItems_;
};

#endif // BPTREECHILDNODE_H
