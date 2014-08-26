#ifndef RPTREECHILDNODE_H
#define RPTREECHILDNODE_H



// ***************************************************************************
//   RPTreeChildNode.h (c) 2014
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

#include "RPTreeNode.h"
#include "RPTreeChildNodeItem.h"
#include <vector>
#include <string>

class RPTreeChildNode : public RPTreeNode
{
    public:
        RPTreeChildNode();
        virtual ~RPTreeChildNode();

        RPChromosomeRegion* getChromosomeBounds();
        int32_t compareRegions(RPChromosomeRegion* chromosomeRegion);
        bool isLeaf();
        int32_t getItemCount();
        RPTreeNodeItem* getItem(int32_t index);
        bool insertItem(RPTreeNodeItem* item);
        bool deleteItem(int32_t index);
    protected:
    private:

        RPChromosomeRegion* chromosomeBounds_;  // chromosome bounds for entire node
        std::vector<RPTreeChildNodeItem* > childItems_; // array for child items

};

#endif // RPTREECHILDNODE_H
