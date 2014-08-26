#ifndef BPTREENODE_H_INCLUDED
#define BPTREENODE_H_INCLUDED

// ***************************************************************************
//   BPTreeNode.h (c) 2014
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


#include "stdint.h"
#include <string>
/**< Interface */

class BPTreeNodeItem;

class BPTreeNode {

    /*
    *   Method returns the node index in the B+ tree organization.
    *
    *   Returns:
    *       node index in B+ tree
    * */
public:

    virtual ~BPTreeNode() {}

    virtual int64_t getNodeIndex() = 0;

    /*
    *   Method identifies the node as a leaf node or a child (non-leaf) node.
    *
    *   Returns:
    *       true, if leaf node; false if child node
    * */
    virtual bool isLeaf()= 0;

    /*
    *   Method inserts the node item appropriate to the item's key value.
    *
    *   Returns:
    *       Node item inserted successfully.
    * */
    virtual bool insertItem(BPTreeNodeItem* item)= 0;

    /*
    *   Method deletes the node item appropriate to the item's index.
    *
    *   Returns:
    *       Node item deleted successfully.
    * */
    virtual bool deleteItem(int index)= 0;

    /*
    *   Method returns the number of items assigned to the node.
    *
    *   Returns:
    *       Count of node items contained in the node
    * */
    virtual int32_t getItemCount()= 0;

    /*
    *   Method returns the indexed node item.
    *
    *   Returns:
    *       Indexed node item.
    * */
    virtual BPTreeNodeItem* getItem(int index)= 0;

    /*
    *   Method returns the lowest chromosome name key belonging to the node.
    *
    *   Returns:
    *       Lowest contig/chromosome name key; or null for no node items.
    * */
    virtual std::string getLowestChromKey()= 0;

    /*
    *   Method returns the highest chromosome name key belonging to the node.
    *
    *   Returns:
    *       Highest contig/chromosome name key; or null for no node items.
    * */
    virtual std::string getHighestChromKey()= 0;

    /*
    *   Method returns the lowest chromosome ID belonging to the node.
    *
    *   Returns:
    *       Lowest contig/chromosome ID; or -1 for no node items.
    * */
    virtual int32_t getLowestChromID()= 0;

    /*
    *   Method returns the highest chromosome ID  belonging to the node.
    *
    *   Returns:
    *       Highest contig/chromosome ID; or -1 for no node items.
    * */
    virtual int32_t getHighestChromID()= 0;

};


#endif // BPTREE_H_INCLUDED
