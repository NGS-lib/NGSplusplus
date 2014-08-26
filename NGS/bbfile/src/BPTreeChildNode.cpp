
// ***************************************************************************
//   BPTreeChildNOde.cpp (c) 2014
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
#include "BPTreeChildNode.h"

BPTreeChildNode::BPTreeChildNode()
{
    //ctor
    isLeafNode_=false;
}

BPTreeChildNode::~BPTreeChildNode()
{
    //dtor
}


 /*
    *   Constructor for the B+ tree child (non-leaf) node.
    *
    *   Parameters:
    *       nodeIndex - index assigned to the node
    *       parent - parent node (object)
    *
    *   Note: Inserted child items contain child/leaf nodes assigned.
    * */
    BPTreeChildNode::BPTreeChildNode(int64_t nodeIndex){

        nodeIndex_ = nodeIndex;
        isLeafNode_ =false;
    }

    /*
    *   Method inserts the node item appropriate to the item's key value.
    *
    *   Returns:
    *       Node item inserted successfully.
    * */
    bool BPTreeChildNode::insertItem(BPTreeNodeItem* item){

        // Quick implementation: assumes all keys are inserted in rank order
        // todo: verify if need to compare key and insert at rank location
        childItems_.push_back((BPTreeChildNodeItem* )item );

        BPTreeNode* childNode = ((BPTreeChildNodeItem* )item)->getChildNode();

        // Note: assumes rank order insertions
        if(childItems_.size() == 1 ){
            lowestChromKey_ = childNode->getLowestChromKey();
            lowestChromID_ = childNode->getLowestChromID();
        }
        else {
            highestChromKey_ = childNode->getHighestChromKey();
            highestChromID_ = childNode->getHighestChromID();
        }

        return true;    // success
    }

    /*
    *   Method deletes the node item appropriate to the item's index.
    *
    *   Returns:
    *       Node item deleted successfully.
    * */
    bool BPTreeChildNode::deleteItem(int index){

        BPTreeChildNodeItem* toErase= childItems_.at(index);
        childItems_.erase(childItems_.begin()+index);
        delete(toErase);
        return true;
    }


    /*
    *   Method returns the indexed node item.
    *
    *   Returns:
    *       node index in B+ tree
    * */
    BPTreeNodeItem* BPTreeChildNode::getItem(int32_t index){
        int itemCount = getItemCount();

        if(index >= itemCount)
            return NULL;

        return childItems_.at(index);
    }

     int32_t BPTreeChildNode::getItemCount() {
        return childItems_.size();
    }


    /*
    *   Method returns the lowest chromosome key value belonging to the node.
    *
    *   Returns:
    *       Lowest contig/chromosome name key value; or null if no node items
    * */
    std::string BPTreeChildNode::getLowestChromKey(){
        if(childItems_.size() > 0)
            return lowestChromKey_;
        else
            return "";
    }

    /*
    *   Method returns the highest chromosome key value belonging to the node.
    *
    *   Returns:
    *       Highest contig/chromosome name key value; or null if no node items
    * */
    std::string BPTreeChildNode::getHighestChromKey(){
        if(childItems_.size() > 0)
            return highestChromKey_;
        else
            return "";
    }

    /*
    *   Method returns the lowest chromosome ID belonging to the node.
    *
    *   Returns:
    *       Lowest key contig/chromosome ID; or -1 if no node items
    * */
    int32_t BPTreeChildNode::getLowestChromID(){
        if(childItems_.size() > 0)
            return lowestChromID_;
        else
            return -1;
    }

    /*
    *   Method returns the highest chromosome ID belonging to the node.
    *
    *   Returns:
    *       Highest key contig/chromosome ID; or -1 if no node items
    * */
     int32_t BPTreeChildNode::getHighestChromID(){
        if(childItems_.size() > 0)
            return highestChromID_;
        else
            return -1;
    }

    /*
    // *********** BPTreeChildNode specific methods *************
    *   Method returns all child items mContained by this child node.
    *
    *   Returns:
    *       List of child items contained by this node
    * */
    std::vector<BPTreeChildNodeItem*>  BPTreeChildNode::getChildItems(){
        return childItems_;
    }
