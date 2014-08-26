// ***************************************************************************
//   ZoomLevelIterator.cpp (c) 2014
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


#include "ZoomLevelIterator.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "BPTree.h"
#include "RPTree.h"
#include "RPChromosomeRegion.h"

ZoomLevelIterator::ZoomLevelIterator()
{
   empty_ = false;
}

ZoomLevelIterator::~ZoomLevelIterator()
{
    //dtor
}


    /**
     * Constructs a zoom level iterator over the specified chromosome region
     * <p/>
     * Parameters:
     * fis - file input stream handle
     * chromIDTree - B+ index tree returns chromId for chromosome name key
     * zoomLevelTree - zoom level R+ chromosome index tree
     * zoomLevel - zoom level represented by the R+ tree
     * selectionRegion_ - chromosome region for selection of Bed feature extraction
     * consists of:
     * startChromID - ID of start chromosome
     * startBase - starting base position for features
     * endChromID - ID of end chromosome
     * endBase - starting base position for features
     * contained - specifies bed features must be contained by region, if true;
     * else return any intersecting region features
     */
     ZoomLevelIterator::ZoomLevelIterator(std::ifstream* fis, BPTree* chromIDTree, RPTree* zoomDataTree_,
                             int32_t zoomLevel, RPChromosomeRegion* selectionRegion_, bool contained)
    {
         empty_ = false;
        // check for valid selection region
        if (selectionRegion_ == NULL)
            throw std::runtime_error("Error: ZoomLevelIterator selection region is null\n");

        this->fis_ = fis;
        this->chromIDTree_ = chromIDTree;
        this->zoomDataTree_ = zoomDataTree_;
        this->zoomLevel_ = zoomLevel;
        this->selectionRegion_ = selectionRegion_;
        isContained_ = contained;

        // set up hit list and read in the first data block
        int32_t hitCount = getHitRegion(selectionRegion_, contained);
        if (hitCount == 0) {
            empty_ = true;
        }

        // Ready for next() data extraction
    }

    /*
     *  Method returns status on a "next record" being available.
     *
     *  Return:
     *      true if a "next record" exists; else false.
     *
     *  Note: If "next" method is called for a false condition,
     *      an UnsupportedOperationException will be thrown.
     * */

    bool ZoomLevelIterator::hasNext() {

        if (empty_)
            return false;

        // first check if current data block can be read for next
        if (zoomRecordIndex_ < (int)zoomRecordList_.size())
            return true;

            // need to fetch next data block
        else if (leafItemIndex_ < (int)leafHitList_.size())
            return true;

        else
            return false;
    }

    /**
     * Method returns the current bed feature and advances to the next bed record.
     * <p/>
     * Returns:
     * Bed feature for current BigBed data record.
     * <p/>
     * Note: If "next" method is called when a "next item" does not exist,
     * an UnsupportedOperationException will be thrown.
     */
    ZoomDataRecord* ZoomLevelIterator::next() {

        // Is there a need to fetch next data block?
        if (zoomRecordIndex_ < (int)zoomRecordList_.size()){
                return  (zoomRecordList_.at(zoomRecordIndex_++));
        }

            // attempt to get next leaf item data block
        else {
            int32_t nHits = getHitRegion(selectionRegion_, isContained_);

            if (nHits > 0) {
                // Note: getDataBlock initializes bed feature index to 0

                return  (zoomRecordList_.at(zoomRecordIndex_++));
              //  return (zoomRecordList_.at(zoomRecordIndex_++)); // return 1st Data Block item
            } else {

                return NULL;
                //throw new NoSuchElementException(result);
            }
        }

    }

    // ************ ZoomLevelIterator specific methods *******************
    /*
   *   Method returns the zoom level assigned to the iterator.
   *
   *   Returns:
   *       Number of leaf node hits allowed at a time
   * */

    int32_t ZoomLevelIterator::getZoomLevel() {
        return zoomLevel_;
    }

    /*
    *   Method returns the iterator selection region.
    * */

    RPChromosomeRegion* ZoomLevelIterator::getSelectionRegion() {
        return selectionRegion_;
    }

/*
    *   Method provides the iterator with a new selection region.
    *
    *   Parameters:
    *      selectionRegion_ - chromosome region for selection of Bed feature extraction
    *      consists of:
    *          startChromID - ID of start chromosome
    *          startBase - starting base position for features
    *          endChromID - ID of end chromosome
    *          endBase - starting base position for features
    *      contained - specifies bed features must be contained by region, if true;
    *          else return any intersecting region features
    *
    *   Returns:
    *       number of chromosome regions found in the selection region
    * */

    int32_t ZoomLevelIterator::setSelectionRegion(RPChromosomeRegion* selectionRegion_,
                                  bool contained) {
        this->selectionRegion_ = selectionRegion_;
        isContained_ = contained;

        // set up hit list and first data block read
        leafHitList_.clear();    // Must nullify existing hit list first!
        int32_t hitCount = getHitRegion(selectionRegion_, contained);
        if (hitCount == 0)   // no hits - no point in fetching data
            throw std::runtime_error("No wig data found in the selection region");

        // Ready for next() data extraction

        return hitCount;
    }


    /*
    *   Method returns if bed items must be completely contained in
    *   the selection region.
    *
    *   Returns:
    *       bool indicates items must be contained in selection region if true,
    *       else may intersect the selection region if false
    * */

     bool ZoomLevelIterator::isContained() {
        return isContained_;
    }


    /*
    *   Method finds the chromosome data hit items for the current hit selection region,
    *   and loads first hit data.
    *
    *   Parameters:
    *       hitRegion - selection region for extracting hit items
    *       contained - indicates hit items must contained in selection region if true;
    *       and if false, may intersect selection region
    *   Note: The selection region will be limited to accommodate  mMaxLeafHits; which terminates
    *       selection at the leaf node at which maxLeafHits is reached. Total number of selected
    *       items may exceed maxLeafHits, but only by the number of leaves in the cutoff leaf node.
    *
    *   Returns:
    *       number of R+ chromosome data hits
    * */

    int32_t ZoomLevelIterator::getHitRegion(RPChromosomeRegion* hitRegion, bool contained) {

        int32_t hitCount = 0;

        // check if new hit list is needed
        // Note: getHitList will reset mLeafItemIndex to 0, the beginning of new hit list
        if (leafHitList_.size()==0) {   //|| mLeafItemIndex >= mLeafHitList.size()){
            hitCount = getHitList(hitRegion, contained);
            if (hitCount == 0)
                return 0;   // no hit data found
        } else {
            hitCount = leafHitList_.size() - leafItemIndex_;
            if (hitCount == 0)
                return 0;   // hit list exhausted
        }

        // Perform a block read for starting base of selection region - use first leaf hit
        dataBlockRead_ = getDataBlock(leafItemIndex_++);

        // try next item - probably intersection issue
        // Note: recursive call until a block is valid or hit list exhuasted
        if (!dataBlockRead_)
            hitCount = getHitRegion(hitRegion, contained);

        return hitCount;
    }

    /*
    *   Method finds the R+ chromosome data tree items for the hit region.
    *
    *   Parameters:
    *       hitRegion - selection region for extracting hit items
    *       contained - indicates hit items must contained in selection region if true;
    *       and if false, may intersect selection region
    *
    *   Note: The selection region will be limited to accommodate  mMaxLeafHits; which terminates
    *       selection at the leaf node at which maxLeafHits is reached. Total number of selected
    *       items may exceed maxLeafHits, but only by the number of leaves in the cutoff leaf node.
    *
    *   Returns:
    *       number of R+ chromosome data hits
    * */

    int32_t ZoomLevelIterator::getHitList(RPChromosomeRegion* hitRegion, bool contained) {

        // hit list for hit region; subject to mMaxLeafHits limitation
        leafHitList_ = zoomDataTree_->getChromosomeDataHits(hitRegion, contained);

        // check if any leaf items were selected
        int32_t nHits = leafHitList_.size();
        if (nHits == 0)
            return 0;   // no data hits found
        else
            leafItemIndex_ = 0;    // reset hit item index to start of list

        // find hit bounds
        int32_t startChromID = leafHitList_.at(0)->getChromosomeBounds()->getStartChromID();
        int32_t startBase = leafHitList_.at(0)->getChromosomeBounds()->getStartBase();
        int32_t endChromID = leafHitList_.at(nHits - 1)->getChromosomeBounds()->getEndChromID();
        int32_t endBase = leafHitList_.at(nHits - 1)->getChromosomeBounds()->getEndBase();

        // save hit region definition; not currently used but useful for debug
        this->hitRegion_ = new RPChromosomeRegion(startChromID, startBase, endChromID, endBase);

        return nHits;
    }

/*
    *   Method sets up a decompressed data block of zoom data records for iteration.
    *
    *   Parameters:
    *       leafItemIndex_ - leaf item index in the hit list referencing the data block
    *
    *   Returns:
    *       Successful Zoom data block set up: true or false.
    * */

    bool ZoomLevelIterator::getDataBlock(int32_t leafItemIndex_) {

        // check for valid data block
        if ((int)leafHitList_.size()==0 || leafItemIndex_ >= (int)leafHitList_.size())
            return false;

        // Perform a block read for indexed leaf item
        leafHitItem_ = leafHitList_.at(leafItemIndex_);

        // get the chromosome names associated with the hit region ID's
        int32_t startChromID = leafHitItem_->getChromosomeBounds()->getStartChromID();
        int32_t endChromID = leafHitItem_->getChromosomeBounds()->getEndChromID();
        chromosomeMap_ = chromIDTree_->getChromosomeIDMap(startChromID, endChromID);

       // bool isLowToHigh = zoomDataTree_.isIsLowToHigh();
        int32_t uncompressBufSize = zoomDataTree_->getUncompressBuffSize();

        // decompress leaf item data block for feature extraction
        zoomDataBlock_ = new ZoomDataBlock(zoomLevel_, fis_, leafHitItem_, &chromosomeMap_,
                uncompressBufSize);

        // get data block zoom data record list and set next index to first item
        zoomRecordList_ = zoomDataBlock_->getZoomData(selectionRegion_, isContained_);
        zoomRecordIndex_ = 0;

        // data block items available for iterator
        if (zoomRecordList_.size() > 0)
            return true;
        else
            return false;
    }



