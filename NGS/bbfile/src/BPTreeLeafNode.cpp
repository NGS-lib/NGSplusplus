// ***************************************************************************
//   BPTreeLeafNode.cpp (c) 2014
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

#include "BPTreeLeafNode.h"
#include "BPTreeNodeItem.h"
#include "BPTreeLeafNodeItem.h"
#include "smart_pointer.h"
BPTreeLeafNode::~BPTreeLeafNode()
{
    //dtor
}

 /*
    *   Constructor for the B+ tree leaf (terminal) node.
    *
    *   Parameters:
    *       nodeIndex - index assigned to the node
    *       parent - parent node (object)
    *
    *   Note: Inserted leaf items contain associated name key/chromosome ID.
    * */
    BPTreeLeafNode::BPTreeLeafNode(int64_t nodeIndex){
        nodeIndex_ = nodeIndex;
         isLeafNode_=true;
    }


    /*
    *   Method inserts the node item appropriate to the item's key value.
    *
    *   Returns:
    *       Node item inserted successfully.
    * */
    bool BPTreeLeafNode::insertItem(BPTreeNodeItem* item){

         // Quick implementation: assumes all keys are inserted in rank order
        // todo: verify if need to compare key and insert at rank location
        BPTreeLeafNodeItem *pToAdd =dynamic_cast<BPTreeLeafNodeItem*>(item);
        leafItems_.push_back(pToAdd);

        // Note: assumes rank order insertions
        if(leafItems_.size() == 1 ){
            lowestChromKey_ = item->getChromKey();
            lowestChromID_ =  dynamic_cast<BPTreeLeafNodeItem*>(item)->getChromID();
        }
        else {
           highestChromKey_ = item->getChromKey();
           highestChromID_ =  dynamic_cast<BPTreeLeafNodeItem*>(item)->getChromID();
        }

        return true;
    }

    /*
    *   Method deletes the node item appropriate to the item's index.
    *
    *   Returns:
    *       Node item deleted successfully.
    * */
    bool BPTreeLeafNode::deleteItem(int32_t index){

        // unacceptable index
        if(index < 0 || index >= getItemCount())
           return false;


        BPTreeLeafNodeItem* toErase= leafItems_.at(index);
        leafItems_.erase(leafItems_.begin()+index);
        delete(toErase);
        return true;  // success
    }

    /*
    *   Method returns the number of items assigned to the node.
    *
    *   Returns:
    *       Count of node items contained in the node
    * */
    int32_t BPTreeLeafNode::getItemCount() {
        return leafItems_.size();
    }

    /*
    *   Method returns the indexed node item.
    *
    *   Returns:
    *       Indexed node item.
    * */
     BPTreeNodeItem* BPTreeLeafNode::getItem(int32_t index){
        if(getItemCount() > 0 && index < getItemCount())
            return leafItems_.at(index);

        return NULL;
    }

    // ************** BPTreeLeafNode specific methods ***********
    /*
    *   Method returns all leaf items mContained by this leaf node.
    *
    *   Returns:
    *       List of leaf items contained by this node
    * */
     std::vector<BPTreeLeafNodeItem*> BPTreeLeafNode::getLeafItems() {
        return leafItems_;
    }
