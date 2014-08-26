// ***************************************************************************
//   RPChromosomeRegion.cpp (c) 2014
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

#include "RPChromosomeRegion.h"

    RPChromosomeRegion::RPChromosomeRegion(){}

    RPChromosomeRegion::~RPChromosomeRegion()
    {
        //dtor
    }


    RPChromosomeRegion::RPChromosomeRegion(uint32_t startChromID, uint32_t startBase,
                              uint32_t endChromID, uint32_t endBase) :startChromID_(startChromID),startBase_(startBase), endChromID_(endChromID),endBase_(endBase)
    {


    }
    /*
   *   Construct region from an existing region.
   * */
    RPChromosomeRegion::RPChromosomeRegion(const RPChromosomeRegion* region) {

        startChromID_ = region->startChromID_;
        startBase_ = region->startBase_;
        endChromID_ = region->endChromID_;
        endBase_ = region->endBase_;
    }


  /*
    *   Comparator for mChromosome bounds is used to find relevant intervals and
    *   rank placement of node items. Returned value indicates relative
    *   positioning to supplied chromosome test region , and expands on normal
    *   comparator by indicating partial overlap in the extremes.
    *
    *   Returns:
    *       - 2 indicates that this region is completely disjoint below the test region
    *       -1 indicates this region intersects the test region from below
    *       0 indicates that this region is inclusive to the test region
    *       1 indicates this region intersects the test region from above
    *       2 indicates that this region is completely disjoint above the test region
    *
    *   Note: additional tests can be applied to determine intersection from above
    *   or below the test region and disjoint above or below the test region cases.
    * */

    uint32_t RPChromosomeRegion::compareRegions(RPChromosomeRegion* testRegion) {

        // test if this region is contained by (i.e. subset of) testRegion region
        if (this->containedIn(testRegion))
            return 0;

            // test if  testRegion region is disjoint from above or below
        else if (this->disjointBelow(testRegion))
            return -2;
        else if (this->disjointAbove(testRegion))
            return 2;

            // Otherwise this region must intersect
        else if (this->intersectsBelow(testRegion))
            return -1;
        else if (this->intersectsAbove(testRegion))
            return 1;

        // unexpected condition is unknown
        return 3;
    }

    /*
    *   Method checks if test region matches this region
    *
    *   Parameters:
    *       testRegion - chromosome selection region
    *
    *   Returns:
    *       This region equals the test region: true or false
    * */

    bool RPChromosomeRegion::equals(RPChromosomeRegion* testRegion) {

        if (startChromID_ == testRegion->startChromID_ && startBase_ == testRegion->startBase_ &&
                endChromID_ == testRegion->endChromID_ && endBase_ == testRegion->endBase_)
            return true;
        else
            return false;
    }

    /*
    *   Method checks if test region contains this region;
    *   (i.e this region is subset oftest region).
    *
    *   Parameters:
    *       testRegion - chromosome selection region
    *
    *   Returns:
    *       This region is contained in the test region: true or false
    * */

    bool RPChromosomeRegion::containedIn(RPChromosomeRegion* testRegion) {

        if (startChromID_ > testRegion->startChromID_ ||
                (startChromID_ == testRegion->startChromID_ && startBase_ >= testRegion->startBase_)) {
            if (endChromID_ < testRegion->endChromID_ ||
                    (endChromID_ == testRegion->endChromID_ && endBase_ <= testRegion->endBase_))
                return true;
            else
                return false;
        } else
            return false;
    }

    /*
    *   Method checks if this region intersects test region from below
    *
    *   Note: To be true, this region must have some part outside the test region
    *
    *   Parameters:
    *       testRegion - chromosome selection region
    *
    *   Returns:
    *       This region intersects the test region from below: true or false
    * */

    bool RPChromosomeRegion::intersectsBelow(RPChromosomeRegion* testRegion) {

        // Only need to test if some part of this region is below and some within test region.
        if (startChromID_ < testRegion->startChromID_ ||
                (startChromID_ == testRegion->startChromID_ && startBase_ < testRegion->startBase_)) {
            if (endChromID_ > testRegion->startChromID_ ||
                    (endChromID_ == testRegion->startChromID_ && endBase_ > testRegion->startBase_))
                return true;
            else
                return false;
        } else
            return false;
    }

    /*
    *   Method checks if this region intersects test region from above.
    *
    *   Note: To be true, this region must have some part outside the test region
    *
    *   Parameters:
    *       testRegion - chromosome selection region
    *
    *   Returns:
    *       This region intersects the test region from above: true or false
    * */

    bool RPChromosomeRegion::intersectsAbove(RPChromosomeRegion* testRegion) {

        // Only need to test if some part of this region is above and some within test region.
        if (endChromID_ > testRegion->endChromID_ ||
                (endChromID_ == testRegion->endChromID_ && endBase_ > testRegion->endBase_)) {
            if (startChromID_ < testRegion->endChromID_ ||
                    (startChromID_ == testRegion->endChromID_ && startBase_ < testRegion->endBase_))
                return true;
            else
                return false;
        } else
            return false;
    }

    /*
    *   Method checks if this region is completely below test region.
    *
    *   Parameters:
    *       testRegion - chromosome selection region
    *
    *   Returns:
    *       This region is disjoint below the test region: true or false
    * */

    bool RPChromosomeRegion::disjointBelow(RPChromosomeRegion* testRegion) {

        if ( (endChromID_ < testRegion->startChromID_) ||
                 ( (endChromID_ == testRegion->startChromID_) && (endBase_ <= testRegion->startBase_) ) )
            return true;
        else
            return false;
    }

    /*
    *   Method checks if this region region is completely above test region.
    *
    *   Parameters:
    *       testRegion - chromosome selection region
    *
    *   Returns:
    *       This region is disjoint above the test region: true or false
    * */

    bool RPChromosomeRegion::disjointAbove(RPChromosomeRegion* testRegion) {

        if ( (startChromID_ > testRegion->endChromID_) ||
                ( startChromID_ == testRegion->endChromID_ && startBase_ >= testRegion->endBase_) )
            return true;
        else
            return false;
    }

    /*
    *   Method computes the extremes between this region and the test region
    *
    *   Parameters:
    *       testRegion - chromosome region to compare against this region
    *
    *   Returns:
    *       new chromosome region of extremes
    * */

    RPChromosomeRegion* RPChromosomeRegion::getExtremes(RPChromosomeRegion* testRegion) {

       RPChromosomeRegion* newRegion = new RPChromosomeRegion(this->getStartChromID(),this->getStartBase(),this->getEndChromID(),this->getEndBase());

       // RPChromosomeRegion* newRegion = new RPChromosomeRegion(*this);

        // update node bounds

        if (testRegion->startChromID_ < newRegion->startChromID_ ||
                (testRegion->startChromID_ == newRegion->startChromID_ &&
                        testRegion->startBase_ < newRegion->startBase_)) {
            newRegion->startChromID_ = testRegion->startChromID_;
            newRegion->startBase_ = testRegion->startBase_;
        }

        if (testRegion->endChromID_ > newRegion->endChromID_ ||
                (testRegion->endChromID_ == newRegion->endChromID_ &&
                        testRegion->endBase_ > newRegion->endBase_)) {
            newRegion->endChromID_ = testRegion->endChromID_;
            newRegion->endBase_ = testRegion->endBase_;
        }

        return newRegion;
    }

    void RPChromosomeRegion::expand(RPChromosomeRegion* testRegion) {
        if (testRegion->startChromID_ < startChromID_ ||
                (testRegion->startChromID_ == startChromID_ &&
                        testRegion->startBase_ < startBase_)) {
            startChromID_ = testRegion->startChromID_;
            startBase_ = testRegion->startBase_;
        }

        if (testRegion->endChromID_ > endChromID_ ||
                (testRegion->endChromID_ == endChromID_ &&
                        testRegion->endBase_ > endBase_)) {
            endChromID_ = testRegion->endChromID_;
            endBase_ = testRegion->endBase_;
        }


    }
