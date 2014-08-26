// ***************************************************************************
//   RpTreeLeafNode.cpp (c) 2014
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

#include "RPTreeLeafNode.h"
#include <cstddef>
RPTreeLeafNode::~RPTreeLeafNode()
{
    if  ( chromosomeBounds_ != NULL){
        delete(chromosomeBounds_);
    }
}


   RPTreeLeafNode::RPTreeLeafNode(){
        // init with null bounds
        chromosomeBounds_ = new RPChromosomeRegion();
    }

    bool RPTreeLeafNode::isLeaf() {
        return true;
    }

    RPChromosomeRegion* RPTreeLeafNode::getChromosomeBounds(){
         return chromosomeBounds_;
    }

    int32_t RPTreeLeafNode::compareRegions(RPChromosomeRegion* chromosomeRegion){

        int32_t value = chromosomeBounds_->compareRegions(chromosomeRegion);
        return value;
    }

    int32_t RPTreeLeafNode::getItemCount() {
        return leafItems_.size();
    }



    RPTreeNodeItem* RPTreeLeafNode::getItem(int32_t index){

       if(index < 0 || index >= (int)leafItems_.size())
            return NULL;
       else
            return leafItems_.at(index);
    }

    bool RPTreeLeafNode::insertItem(RPTreeNodeItem* item){

         RPTreeLeafNodeItem* newItem =  (RPTreeLeafNodeItem*)item;

        // Note: assumes all keys are inserted in rank order
        leafItems_.push_back(newItem);

        // todo: compare region and insert at appropriate indexed rank location
        //   leafHitItem.add( index, (RPTreeLeafNodeItem)item );

        // update leaf node chromosome bounds - use extremes
        // Update node bounds or start node chromosome bounds with first entry
       if(chromosomeBounds_ == NULL)
            chromosomeBounds_ = new RPChromosomeRegion(newItem->getChromosomeBounds());
       else
            chromosomeBounds_ = chromosomeBounds_->getExtremes(newItem->getChromosomeBounds());

        // successful insert
         return true;
    }

     bool RPTreeLeafNode::deleteItem(int32_t index){

        int32_t itemCount = getItemCount();

        // unacceptable index  - reject
        if(index < 0 || index >= itemCount)
            return false;

        // delete indexed entry
        RPTreeLeafNodeItem* toErase= leafItems_.at(index);
        leafItems_.erase(leafItems_.begin()+index);
        delete(toErase);

        // successful delete
        return true;
    }

