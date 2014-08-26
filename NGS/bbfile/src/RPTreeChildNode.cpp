// ***************************************************************************
//   RPTreeChildNode.cpp (c) 2014
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

#include "RPTreeChildNode.h"

RPTreeChildNode::~RPTreeChildNode()
{
    //dtor
}


    RPTreeChildNode::RPTreeChildNode() {
        chromosomeBounds_=NULL;
        // Note: Chromosome bounds are null until a valid region is specified
    }

    // *** BPTreeNode interface implementation ***


    RPChromosomeRegion* RPTreeChildNode::getChromosomeBounds() {
        return chromosomeBounds_;
    }

     int32_t RPTreeChildNode::compareRegions(RPChromosomeRegion* chromosomeRegion) {

        // test leaf item bounds for hit
        int32_t value = chromosomeBounds_->compareRegions(chromosomeRegion);
        return value;
    }

    bool RPTreeChildNode::isLeaf() {
        return false;
    }

    int32_t RPTreeChildNode::getItemCount() {
        return childItems_.size();
    }

    RPTreeNodeItem* RPTreeChildNode::getItem(int32_t index) {

        if (index < 0 || index >= (int)childItems_.size())
            return NULL;
        else {
            RPTreeChildNodeItem* item = childItems_.at(index);
            return (RPTreeNodeItem*) item;
        }
    }

    bool RPTreeChildNode::insertItem(RPTreeNodeItem* item) {

        RPTreeChildNodeItem* newItem = (RPTreeChildNodeItem*) item;

        // Quick implementation: assumes all keys are inserted in rank order
        // todo: or compare key and insert at rank location
        childItems_.push_back(newItem);

        // Update node bounds or start node chromosome bounds with first entry
        if (chromosomeBounds_ == NULL){
            chromosomeBounds_ = new RPChromosomeRegion(newItem->getChromosomeBounds());
            }
        else
        {
              chromosomeBounds_ = chromosomeBounds_->getExtremes(newItem->getChromosomeBounds());
        }


        // success
        return true;
    }




    bool RPTreeChildNode::deleteItem(int32_t index) {

        int32_t itemCount = getItemCount();

        // unacceptable index  - reject
        if (index < 0 || index >= itemCount)
            return false;

        // delete indexed entry
        RPTreeChildNodeItem* toErase= childItems_.at(index);
        childItems_.erase(childItems_.begin()+index);
        delete(toErase);
        // successful delete
        return true;
    }

