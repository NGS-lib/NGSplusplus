#ifndef BIGWIGSECTION_H
#define BIGWIGSECTION_H

// ***************************************************************************
//   BigWigSection.h (c) 2014
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


#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "stdint.h"
#include "RPTreeLeafNodeItem.h"
#include "BigWigSectionHeader.h"
class BigWigSection
{
    public:
        BigWigSection();
        virtual ~BigWigSection();
        BigWigSection(std::vector<char>& sectionBuffer, std::map<uint32_t,std::string>* chromosomeMap,
                         std::vector<RPTreeLeafNodeItem*>::iterator leafHitItem);
        bool isValidSectionType();
        int32_t getItemCount();
        BigWigSectionHeader* getSectionHeader();
        int32_t getSectionDataSize();
        int32_t getSectionData(RPChromosomeRegion* selectionRegion, bool contained,
                              std::vector<WigItem*>& wigItemList);



    protected:
    private:
        std::stringstream myStream_;
      //std::stringstream myStream;  std::ifstream* dis_;

        std::vector<RPTreeLeafNodeItem*>::iterator leafHitItem_;
        int32_t sectionDataSize_;
        std::map<uint32_t , std::string>* chromosomeMap_;
        BigWigSectionHeader* wigSectionHeader_;
        std::vector<char> sectionBuffer_;
};

#endif // BIGWIGSECTION_H
