#ifndef ZOOMLEVELITERATOR_H
#define ZOOMLEVELITERATOR_H


// ***************************************************************************
//   ZoomLevelIterator.h (c) 2014
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




#include <map>
#include <vector>
#include <string>
#include "ZoomDataBlock.h"
#include "ZoomDataRecord.h"
#include "BPTree.h"
#include "RPTree.h"


class ZoomLevelIterator
{
    public:
        ZoomLevelIterator();
        virtual ~ZoomLevelIterator();

      ZoomLevelIterator(std::ifstream* fis, BPTree* chromIDTree, RPTree* zoomDataTree,
                             int32_t zoomLevel, RPChromosomeRegion* selectionRegion, bool contained);
     bool hasNext();
     ZoomDataRecord* next();
     int32_t getZoomLevel();
     RPChromosomeRegion* getSelectionRegion();
     int32_t setSelectionRegion(RPChromosomeRegion* selectionRegion,
                                  bool contained);
    bool isContained();

    int32_t getHitRegion(RPChromosomeRegion* hitRegion, bool contained);
    int32_t getHitList(RPChromosomeRegion* hitRegion, bool contained);
    bool getDataBlock(int32_t leafItemIndex);

    protected:
    private:

    bool empty_;

    // zoom level for zoom data
    int32_t zoomLevel_;

    //specification of chromosome selection region
    RPChromosomeRegion* selectionRegion_;  // selection region for iterator
    bool isContained_; // if true, features must be fully contained by extraction region
    RPChromosomeRegion* hitRegion_;  // hit selection region for iterator

    // File access variables for reading zoom level data block
    std::ifstream* fis_;  // file input stream handle
    BPTree* chromIDTree_;    // B+ chromosome index tree
    RPTree* zoomDataTree_;  // R+ zoom data locations tree

    // chromosome region extraction items
    std::vector<RPTreeLeafNodeItem*> leafHitList_; // array of leaf hits for selection region items
    std::map<uint32_t, std::string> chromosomeMap_;  // map of chromosome ID's and corresponding names
    int32_t leafItemIndex_;   // index of current leaf item being processed from leaf hit list
    RPTreeLeafNodeItem* leafHitItem_;   // leaf item being processed by next

    // current zoom level block being processed
    ZoomDataBlock* zoomDataBlock_;  // holds data block of zoom level records decompressed
    bool dataBlockRead_;  // flag indicates successful read of data block for current leaf item
    std::vector<ZoomDataRecord*> zoomRecordList_; // array of selected zoom data records
    int32_t zoomRecordIndex_;    // index of next zoom data record from the list

};

    class EmptyIterator:public ZoomLevelIterator
    {

      //  static EmptyIterator theInstance = new EmptyIterator();

        bool hasNext() {
            return false;
        }
    };

#endif // ZOOMLEVELITERATOR_H
