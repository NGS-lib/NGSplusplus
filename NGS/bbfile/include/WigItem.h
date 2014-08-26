#ifndef WIGITEM_H_INCLUDED
#define WIGITEM_H_INCLUDED


// ***************************************************************************
//   WigItem.h (c) 2014
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


#include <string>
#include <vector>
#include "stdint.h"

class WigItem {

private:
     int32_t itemIndex_;         // wig section item index number
     std::string chromosome_;     // mChromosome name
     int32_t startBase_;         // mStartBase base position for feature
     int32_t endBase_;           // mEndBase base position for feature
     float wigValue_;        // wig value

    public:

     WigItem(int32_t itemIndex, std::string chromosome, int32_t startBase, int32_t endBase, float wigValue){

        this->itemIndex_ = itemIndex;
        this->chromosome_ = chromosome;
        this->startBase_ = startBase;
        this->endBase_ = endBase;
        this->wigValue_ = wigValue;
    }

     int32_t getItemNumber(){
        return itemIndex_;
    }

     std::string getChromosome() {
        return chromosome_;
    }

     int32_t getStartBase() {
        return startBase_;
    }

     int32_t getEndBase() {
        return endBase_;
    }

     float getWigValue() {
        return wigValue_;
    }

};


#endif // WIGITEM_H_INCLUDED
