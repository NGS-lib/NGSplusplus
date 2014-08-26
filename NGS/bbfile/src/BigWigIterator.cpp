// ***************************************************************************
//   BigWigIterator.cpp (c) 2014
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



//TODO This entire structure, including the other iterators, needs to be redone
/*
    Currently, the data is stored in the iterator, with the iterator handling
    both the data handling and data access functionalities.

    We need to create an intermediate structure between BBFileReader and the iterator ( DataRange ) that will
    handle storing and buffering the required ranged. The iterator will then only handle access duties.

*/

#include "BigWigIterator.h"
#include <stdexcept>

BigWigIterator::BigWigIterator():empty_(false)
{
   // this->wigDataBlock_==NULL;
//    this->hitRegion_=NULL;
}

BigWigIterator::~BigWigIterator()
{
    //dtor
  //  this->hitRegion_!=NULL;
  //  delete(this->hitRegion_);

   // if (this->wigDataBlock_!=NULL)
  //      delete (this->wigDataBlock_);
}


BigWigIterator::BigWigIterator(std::ifstream* fis, BPTree* chromIDTree, RPTree* chromDataTree,
                               RPChromosomeRegion* selectionRegion, bool contained):
    empty_(false),selectionRegion_(selectionRegion),fis_(fis),isContained_(contained),chromIDTree_(chromIDTree),chromDataTree_(chromDataTree)
{
    try
    {
        // check for valid selection region
        if (selectionRegion == NULL)
            throw std::runtime_error("Error: BigWigIterator selection region is null\n");

        isContained_ = contained;

        // set up hit list and first data block read
        int32_t hitCount = loadNextLeaf(selectionRegion, contained);
        if (hitCount == 0)
        {
            empty_ = true;
        }
    }
    catch(...)
    {
        throw;
    }

}


/**< True if the iterator is pointing to the last element in the Range */
bool BigWigIterator::isEnd()
{
    //If moving on current data block, obviously not finished
    if (wigcurrent_ != wigend_)
    {
        return false;
    }
    //Any other data blocks?
    if ( leafcurrent_!=leafend_)
        return true;

    return true;

}



/** \brief Access next element. Will change once we decouple iterator and data handling
 *
 * \return const BigWigIterator&
 *
 */
const BigWigIterator& BigWigIterator::operator++()
{

    /**< Increment vector block */
    ++wigcurrent_;
    if (wigcurrent_ != wigend_)
    {
        return *this;
    }

    /**< Else increment leaf */
    if (leafcurrent_!=leafend_)
    {
        //TODO recode this
        /**< This sets our iterators correctly..and is horrible horrible horrible way to to do */
        int32_t nHits = loadNextLeaf(selectionRegion_, isContained_);
        if (nHits > 0)
        {
            return *this;

        }
        else
        {
            leafcurrent_=leafend_;
            return *this;
        }

    }
    return *this;


}

BigWigIterator BigWigIterator::operator++(int)
{

    BigWigIterator old(*this);
    ++*this;
    return old;
}

// ************ BigBedIterator specific methods *******************
/*
*   Method returns the iterator selection region.
* */

RPChromosomeRegion* BigWigIterator::getSelectionRegion()
{
    return selectionRegion_;
}

/*
    *   Method provides the iterator with a new selection region.
    *
    *   Parameters:
    *      selectionRegion - chromosome region for selection of Bed feature extraction
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

int32_t BigWigIterator::setSelectionRegion(RPChromosomeRegion* selectionRegion,
        bool contained)
{
    this->selectionRegion_ = selectionRegion;
    isContained_ = contained;

    // set up hit list and first data block read
    //TODO  Surely this leaks

    leafHitList_.clear();// Must nullify existing hit list first!
    int32_t hitCount = loadNextLeaf(selectionRegion, contained);
    if (hitCount == 0)   // no hits - no point in fetching data
        throw new std::runtime_error("No wig data found in the selection region");

    // Ready for next() data extraction
    return hitCount;
}

/*
*   Method returns if bed items must be completely contained in
*   the selection region.
*
*   Returns:
*       Boolean indicates items must be contained in selection region if true,
*       else may intersect the selection region if false
* */

bool BigWigIterator::isContained()
{
    return isContained_;
}

/*
*   Method returns the B+ chromosome index tree used for identifying
*   chromosome ID's used to specify R+ chromosome data locations.
*
*   Returns:
*       B+ chromosome index tree
* */

BPTree* BigWigIterator::getChromosomeIDTree()
{
    return chromIDTree_;
}

/*
*   Method returns the R+ chromosome data tree used for identifying
*   chromosome data locations for the selection region.
*
*   Returns:
*       R+ chromosome data locations tree
* */

RPTree* BigWigIterator::getChromosomeDataTree()
{
    return chromDataTree_;
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
//getHitRegion
int32_t BigWigIterator::loadNextLeaf(RPChromosomeRegion* hitRegion, bool contained)
{
    try
    {
        int32_t hitCount = 0;

        // check if new hit list is needed
        if (leafHitList_.size()==0 )
        {
            hitCount = filterLeafVector(hitRegion, contained);
            if (hitCount == 0)
                return 0;   // no hit data found
        }
        else
        {
            hitCount = std::distance(leafcurrent_,leafend_) ;
            if (hitCount == 0)
                return 0;   // hit list exhausted
        }

        // Extract current leaf and return number of elements in leaf
        dataBlockRead_ = getDataBlock(leafcurrent_);
        ++leafcurrent_;

        //TODO check why this is needed. We should be able to pre-filter a leaf?
        // try next item - probably intersection issue
        // Note: recursive call until a block is valid or hit list exhuasted
        if (!dataBlockRead_)
            hitCount = loadNextLeaf(hitRegion, contained);

        return hitCount;

    }
    catch(...)
    {
        throw;
    }
}

/*
*   Method finds the chromosome data tree hit items for the current hit selection region.
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
* */ //getHitList
int32_t BigWigIterator::filterLeafVector(RPChromosomeRegion* hitRegion, bool contained)
{

    try
    {
        // hit list for hit region; subject to mMaxLeafHits limitation
        leafHitList_ = chromDataTree_->getChromosomeDataHits(hitRegion, contained);

        // check if any leaf items were selected

        leafcurrent_ =  leafHitList_.begin();
        leafend_ = leafHitList_.end();


        int32_t nHits = leafHitList_.size();

        if (nHits == 0)
            return 0;

        // find hit bounds from first and last hit items
    //    int32_t startChromID = leafHitList_.at(0)->getChromosomeBounds()->getStartChromID();
    //    int32_t startBase = leafHitList_.at(0)->getChromosomeBounds()->getStartBase();
    //    int32_t endChromID = leafHitList_.at(nHits - 1)->getChromosomeBounds()->getEndChromID();
     //   int32_t endBase = leafHitList_.at(nHits - 1)->getChromosomeBounds()->getEndBase();

        // save hit region; not currently used but useful for debug
        //this->hitRegion_ = RPChromosomeRegion(startChromID, startBase, endChromID, endBase);

        return nHits;
    }
    catch(...)
    {
        throw;
    }
}

/*
*   Method sets up a decompressed data block of big bed features for iteration.
*
*   Parameters:
*       leafIter iterator pointing to leaf who's block we are doing to decompress
*
*   Returns:
*       Success of operation
* */
bool BigWigIterator::getDataBlock(std::vector<RPTreeLeafNodeItem*>::iterator leafIter)
{

    try
    {
        //leafHitItem_=leafIter;
       // RPTreeLeafNodeItem* test =*leafIter;
        int32_t startChromID = (*leafIter)->getChromosomeBounds()->getStartChromID();
        int32_t endChromID = (*leafIter)->getChromosomeBounds()->getEndChromID();
        chromosomeMap_ = chromIDTree_->getChromosomeIDMap(startChromID, endChromID);

        int32_t uncompressBufSize = chromDataTree_->getUncompressBuffSize();
        // decompress leaf item data block for feature extraction
        //TODO verify this delete
        //Clear Data block? I Guess?
//        if (wigDataBlock_!=NULL)
//            delete (wigDataBlock_);

        wigDataBlock_ = BigWigDataBlock(fis_, leafIter, uncompressBufSize);

        // get section Wig item list and set next index to first item
        wigItemList_ = wigDataBlock_.getWigData(*selectionRegion_,chromosomeMap_, isContained_);
        wigcurrent_ = wigItemList_.begin();
        wigend_ = wigItemList_.end();

        // data block items available for iterator
        if (wigcurrent_ != wigend_)
            return true;
        else
            return false;

    }
    catch(...)
    {
        throw;
    }
}
