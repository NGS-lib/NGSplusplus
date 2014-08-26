#ifndef BEDFEATURE_H
#define BEDFEATURE_H
#include <string>
#include "stdint.h"

// ***************************************************************************
//    BedFeature.h (c) 2014
//    Copyright @ Alexei Nordell-Markovits : Sherbrooke University
//
//    This file is part BWReader library.
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


class BedFeature
{
    public:
        BedFeature() {}
        virtual ~BedFeature() {}
    protected:
    private:

    int32_t itemIndex_;     // data record index
    // BBFile Table I - BigBed data format
     std::string chromosome_;      // mChromosome/contig name
     int32_t startBase_;         // starting base for item
     int32_t endBase_;           // ending base for item
     std::string restOfFields_;    // string containing custom fields

    public:


    BedFeature(int32_t itemIndex, std::string chromosome, int32_t startBase, int32_t endBase, std::string restOfFieldsString){

       this->itemIndex_ = itemIndex;
       this->chromosome_ =  chromosome;
       this->startBase_ =  startBase;
       this->endBase_ = endBase;
    /**< TODOOOOOO! */
    //todo fix this?
     //  restOfFields_ = ( restOfFieldsString == NULL ? NULL : restOfFieldsString.split("\t") );
   }

   // returns the data record index
   int32_t getItemIndex() {
       return itemIndex_;
   }

   // returns the mChromosome ID (0, 1, etc.)
   std::string getChromosome() {
       return chromosome_;
   }

   // returns the mChromosome mStartBase base position
   int32_t getStartBase(){
       return startBase_;
   }

   // returns the mChromosome mEndBase base position
   int32_t getEndBase() {
       return endBase_;
   }

    std::string getRestOfFields(){
        return restOfFields_;
    }

};




#endif // BEDFEATURE_H
