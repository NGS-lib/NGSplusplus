#ifndef BIGWIGITERATOR_H
#define BIGWIGITERATOR_H


// ***************************************************************************
//   BigWigIterator.h (c) 2014
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


#include "WigItem.h"
#include "RPTreeLeafNodeItem.h"
#include "BPTree.h"
#include "RPTree.h"
#include "BigWigDataBlock.h"


#include <vector>
#include <string>
#include <iterator>
#include <map>
#include <fstream>
#include <iostream>

class BigWigIterator:public std::iterator<std::input_iterator_tag, WigItem>
{
    public:
        BigWigIterator();
//        BigWigIterator();
        virtual ~BigWigIterator();

        BigWigIterator( std::ifstream* fis,BPTree* chromIDTree, RPTree* chromDataTree,
                          RPChromosomeRegion* selectionRegion, bool contained);

        bool isEnd();
        BigWigIterator operator++(int);
        const BigWigIterator& operator++();

        //Copy constructor
        BigWigIterator(const BigWigIterator& other) {

             leafHitList_.resize(other.leafHitList_.size());
            for(unsigned int i=0; i<other.leafHitList_.size();i++)
                {
                    leafHitList_.at(i)=other.leafHitList_.at(i);
                }

             wigcurrent_=other.wigcurrent_;
             wigend_=other.wigend_;

            //Set iterator to appropriate distance
             int distance =std::distance(other.leafcurrent_,other.leafend_);
             distance = other.leafHitList_.size()- distance;
             leafcurrent_= (leafHitList_.begin()+ distance);
             leafend_= leafHitList_.end();

             this->empty_= other.empty_;


            selectionRegion_= other.selectionRegion_;
            fis_= other.fis_;
            isContained_= other.isContained_;

            chromIDTree_= other.chromIDTree_;
            chromDataTree_= other.chromDataTree_;


           // leafHitList_= other.leafHitList_;
            chromosomeMap_= other.chromosomeMap_;
            leafItemIndex_= other.leafItemIndex_;

            leafHitItem_= other.leafHitItem_;

            wigDataBlock_= other.wigDataBlock_;
            dataBlockRead_= other.dataBlockRead_;
            wigItemList_= other.wigItemList_;
        }

        void swap(BigWigIterator & other)
        {

             std::swap(wigcurrent_, other.wigcurrent_);
             std::swap(wigend_, other.wigend_);

             std::swap(leafcurrent_, other.leafcurrent_);
             std::swap(leafend_, other.leafend_);

            std::swap(this->empty_, other.empty_);


            std::swap(selectionRegion_, other.selectionRegion_);
            std::swap(fis_, other.fis_);
            std::swap(isContained_, other.isContained_);

            std::swap(chromIDTree_, other.chromIDTree_);
            std::swap(chromDataTree_, other.chromDataTree_);

            std::swap(leafHitList_, other.leafHitList_);
            std::swap(chromosomeMap_, other.chromosomeMap_);
            std::swap(leafItemIndex_, other.leafItemIndex_);

            std::swap(leafHitItem_, other.leafHitItem_);

            std::swap(wigDataBlock_, other.wigDataBlock_);
            std::swap(dataBlockRead_, other.dataBlockRead_);
            std::swap(wigItemList_, other.wigItemList_);

        }

        //Swap trick
        BigWigIterator& operator=(BigWigIterator other)
        {
            this->swap(other);
            return *this;
        }


        RPChromosomeRegion* getSelectionRegion();
        int32_t setSelectionRegion(RPChromosomeRegion* selectionRegion, bool contained);
        int32_t loadNextLeaf(RPChromosomeRegion* hitRegion, bool contained);
        int32_t filterLeafVector(RPChromosomeRegion* hitRegion, bool contained);
        bool getDataBlock(std::vector<RPTreeLeafNodeItem*>::iterator leafIter);


        bool isContained();
        BPTree* getChromosomeIDTree();
        RPTree* getChromosomeDataTree();


         /**< Mandatory overloads for input iterator */

         //Equality operators will not work correctly
      //  friend bool operator ==(const BigWigIterator& cit, const BigWigIterator& cit_cur)
      //  {
      //      return cit_cur.wigcurrent_ == cit.wigcurrent_;
      //  }
      //  friend  bool operator !=(const BigWigIterator& cit, const BigWigIterator& cit_cur)
      //  {
      //      return !(cit_cur == cit);
      //  }
    /**< Complete */
        WigItem operator*() const
        {
            return *wigcurrent_;
        }

        WigItem* operator->() {

            return (&*wigcurrent_);
        }



    protected:
    private:

    std::vector<WigItem>::iterator wigcurrent_;
    std::vector<WigItem>::iterator wigend_;


    std::vector<RPTreeLeafNodeItem*>::iterator leafcurrent_;
    std::vector<RPTreeLeafNodeItem*>::iterator leafend_;


    bool empty_;
    //specification of chromosome selection region
    RPChromosomeRegion* selectionRegion_;  // selection region for iterator
    std::ifstream* fis_;  // file input stream handle
    bool isContained_;     // if true, features must be fully contained by selection region

    // File access variables for reading Bed data block
    BPTree* chromIDTree_;    // B+ chromosome index tree
    RPTree* chromDataTree_;  // R+ chromosome data location tree

    // chromosome region extraction items
    std::vector<RPTreeLeafNodeItem*> leafHitList_; // array of leaf hits for selection region items
    std::map<uint32_t, std::string> chromosomeMap_;  // map of chromosome ID's and corresponding names
    int32_t leafItemIndex_;   // index of current leaf item being processed from leaf hit list

    std::vector<RPTreeLeafNodeItem*>::iterator leafHitItem_;

    // current data block processing members
    BigWigDataBlock wigDataBlock_;  // Wig data block with Wig records decompressed
    bool dataBlockRead_;  // indicates successful read of data block
    std::vector<WigItem> wigItemList_; // array of selected Wig values
//   int32_t wigItemIndex_;      // index of next Wig data item from the list
};

#endif // BIGWIGITERATOR_H
