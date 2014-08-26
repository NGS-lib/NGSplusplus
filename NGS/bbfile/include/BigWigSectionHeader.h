#ifndef BIGWIGSECTIONHEADER_H
#define BIGWIGSECTIONHEADER_H

// ***************************************************************************
//   BigWigSectionHeader.h.h (c) 2014
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
#include "enum_global.h"
class BigWigSectionHeader
{
    public:
        BigWigSectionHeader();
        virtual ~BigWigSectionHeader();

        static const int32_t SECTION_HEADER_SIZE = 24;
        static const int32_t FIXEDSTEP_ITEM_SIZE = 4;
        static const int32_t VARSTEP_ITEM_SIZE = 8;
        static const int32_t BEDGRAPH_ITEM_SIZE = 12;


        BigWigSectionHeader(std::ifstream& bdis);
        BigWigSectionHeader(std::stringstream& bdis);
        uint32_t getChromID() ;
        uint32_t getChromosomeStart();
        uint32_t getChromosomeEnd();
        uint32_t getItemStep();
        uint32_t getItemSpan();
        WigTypeNamespace::WigItemType getItemType();
        bool IsValidType();
        uint16_t getItemCount();
        unsigned char getReserved();
        bool getItemType(char type);



    private:

        uint32_t chromID_;       // Chromosome/contig Numerical ID from BBFile Chromosome B+ tree
        uint32_t chromStart_;    // starting base position
        uint32_t chromEnd_;      // ending base position
        uint32_t itemStep_;      // number of base spaces between fixed items
        uint32_t itemSpan_;      // number of bases in fixed step items
        WigTypeNamespace::WigItemType itemType_; // type of data items: 1 = bedGraph, 2 = varStep, 3 = fixedStep
        unsigned char reserved_;     // reserved; currently = 0
        short itemCount_;   // number of data items in this chromosome section

        bool isValidType_;    // indicates a if a valid Wig item type was read
        std::string itemDescription_; // string representation of item type.

};

#endif // BIGWIGSECTIONHEADER_H
