#ifndef RPTREEE_H
#define RPTREEE_H


// ***************************************************************************
//   RPTree.h (c) 2014
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

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include "stdint.h"
#include "RPChromosomeRegion.h"
#include "RPTreeLeafNodeItem.h"
#include "RPTreeNode.h"
#include "RPTreeHeader.h"
#include "RPTreeChildNodeItem.h"
class RPTree
{
    public:
        RPTree();
        virtual ~RPTree();

        static const  uint32_t RPTREE_NODE_FORMAT_SIZE = 4;       // node format size
        static  const uint32_t RPTREE_NODE_LEAF_ITEM_SIZE = 32;   // leaf item size
        static const  uint32_t RPTREE_NODE_CHILD_ITEM_SIZE = 24;  // child item size


        RPTree(uint32_t order);
        RPTree(std::ifstream& fis, uint64_t fileOffset, uint32_t uncompressBuffSize, bool forceDescend);
        RPChromosomeRegion* getChromosomeRegion(uint32_t startChromID, uint32_t endChromID);
        std::vector<RPChromosomeRegion*> getAllChromosomeRegions();
        RPChromosomeRegion* findChromosomeRegion(RPTreeNode* thisNode, uint32_t startChromID, uint32_t endChromID, RPChromosomeRegion* region);


        static RPTreeNode* readRPTreeNode(std::ifstream* fis, uint64_t fileOffset, bool forceDescend);


        void findAllChromosomeRegions(RPTreeNode* thisNode, std::vector<RPChromosomeRegion* > regionList);
        std::vector<RPTreeLeafNodeItem*> getChromosomeDataHits(RPChromosomeRegion* selectionRegion, bool contained);
        void findChromosomeRegionItems(RPTreeNode* thisNode, RPChromosomeRegion* selectionRegion, std::vector<RPTreeLeafNodeItem*> &leafHitItems);


        uint32_t getUncompressBuffSize() { return uncompressBuffSize_; }
        void setUncompressBuffSize(uint32_t val) { uncompressBuffSize_ = val; }
        uint64_t getRPTreeOffset() { return rpTreeOffset_; }
        void setRPTreeOffset(uint64_t val) { rpTreeOffset_ = val; }
        RPTreeHeader* getRPTreeHeader() { return rpTreeHeader_; }
        void setRPTreeHeader(RPTreeHeader* val) { rpTreeHeader_ = val; }
        RPChromosomeRegion* getChromosomeBounds() { return chromosomeBounds_; }
        void setChromosomeBounds(RPChromosomeRegion* val) { chromosomeBounds_ = val; }
        uint32_t getOrder() { return order_; }
        void setOrder(uint32_t val) { order_ = val; }
        RPTreeNode* getRootNode() { return rootNode_; }
        void setRootNode(RPTreeNode* val) { rootNode_ = val; }
        uint64_t getNodeCount() { return nodeCount_; }
        void setNodeCount(uint64_t val) { nodeCount_ = val; }
        uint64_t getLeafCount() { return leafCount_; }
        void setLeafCount(uint64_t val) { leafCount_ = val; }


        int64_t getItemCount() {
            return rpTreeHeader_->getItemCount();
        }



    protected:
    private:
        uint32_t uncompressBuffSize_;
        uint64_t rpTreeOffset_;
        RPTreeHeader* rpTreeHeader_;
        RPChromosomeRegion* chromosomeBounds_;
        uint32_t order_;
        RPTreeNode* rootNode_;
        uint64_t nodeCount_;
        uint64_t leafCount_;
};

#endif // RPTREEE_H
