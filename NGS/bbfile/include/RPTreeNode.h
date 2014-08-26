#ifndef RPTREENODE_H
#define RPTREENODE_H


// ***************************************************************************
//   RPTreeNode.h (c) 2014
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
#include "RPTreeNodeItem.h"

class RPTreeNode
{
    public:
        RPTreeNode() {}
        virtual ~RPTreeNode() {}
  // Identifies the node as a leaf node or a child (non-leaf) node.
    virtual bool isLeaf()=0;

    // Returns the chromosome bounds belonging to the entire node.
    virtual RPChromosomeRegion* getChromosomeBounds()=0;

     // Note: compareRegions returns the following values:
     //   -2 indicates chromosome region is completely below node region
     //   -1 indicates that chromosome region uint32_tersect node region from below
     //  0 means that chromosome region is inclusive to node region
     //  1 indicates chromosome region uint32_tersects node region from above
     //  2 indicates that this region is completely above that region
    virtual int32_t compareRegions(RPChromosomeRegion* chromosomeRegion)=0;
    // Returns the number of items assigned to the node.
    virtual int32_t getItemCount()=0;
    // Returns the indexed node item.
    virtual RPTreeNodeItem* getItem(int32_t index)=0;
    // Inserts new node item according to bounds rank
    virtual bool insertItem(RPTreeNodeItem* item)=0;
    // Deletes indexed node item
    virtual bool deleteItem(int32_t index)=0;

};

#endif // RPTREENODE_H
