#ifndef RPTREENODEITEM_H_INCLUDED
#define RPTREENODEITEM_H_INCLUDED


// ***************************************************************************
//   RPTreeNodeItem.h (c) 2014
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


#include "RPChromosomeRegion.h"
class RPTreeNodeItem{

    public:
    // returns the chromosome boundary for the item
     virtual RPChromosomeRegion* getChromosomeBounds()=0;

    // Note: compareRegions returns the following values:
     //   -2 indicates chromosome region is completely below node region
     //   -1 indicates that chromosome region intersect node region from below
     //  0 means that chromosome region is inclusive to node region
     //  1 indicates chromosome region intersects node region from above
     //  2 indicates that this region is completely above that region
  //  virtual uint32_t compareRegions(RPChromosomeRegion* chromosomeRegion)=0;
};

#endif // RPTREENODEITEM_H_INCLUDED
