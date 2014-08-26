// ***************************************************************************
//   RPTreeChildNodeItem.cpp (c) 2014
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

#include "RPTreeChildNodeItem.h"
#include "RPChromosomeRegion.h"
#include "RPTreeNodeProxy.h"
#include "RPTree.h"


#include <vector>
#include <string>
RPTreeChildNodeItem::RPTreeChildNodeItem()
{
    //ctor
}

RPTreeChildNodeItem::~RPTreeChildNodeItem()
{
    //dtor
}


  RPTreeChildNodeItem::RPTreeChildNodeItem(uint32_t startChromID, uint32_t startBase,
                               uint32_t endChromID, uint32_t endBase, RPTreeNode* childNode) {


        chromosomeBounds_ = new RPChromosomeRegion(startChromID, startBase, endChromID, endBase);
        this->childNode_ = childNode;
    }


    RPTreeChildNodeItem::RPTreeChildNodeItem(uint32_t startChromID, uint32_t startBase,
                               uint32_t endChromID, uint32_t endBase, uint64_t childDataOffset) {


        chromosomeBounds_ = new RPChromosomeRegion(startChromID, startBase, endChromID, endBase);
    }


    RPChromosomeRegion* RPTreeChildNodeItem::getChromosomeBounds() {
        return chromosomeBounds_;
    }

    RPTreeNode* RPTreeChildNodeItem::getChildNode() {

        RPTreeNodeProxy *p = dynamic_cast<RPTreeNodeProxy*>(childNode_);

        if (p){
             RPTreeNodeProxy* proxy = (RPTreeNodeProxy *) childNode_;
             childNode_ = RPTree::readRPTreeNode( proxy->fis_, proxy->fileOffset_, true);
        }

        return childNode_;
    }

    uint32_t RPTreeChildNodeItem::compareRegions(RPChromosomeRegion* chromosomeRegion) {

        uint32_t value = chromosomeBounds_->compareRegions(chromosomeRegion);
        return value;
    }
