
// ***************************************************************************
//   BBTotalSummaryBlock.cpp (c) 2014
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

#include "BBTotalSummaryBlock.h"
#include "endian_helper.h"
#include <iostream>
#include <fstream>
#include "stdint.h"
#include <stdexcept>

         /*
   *   Constructor for reading in TotalSummaryBlock from BBFile
   *
   *    Parameters:
   *    fis - file input stream handle
   *    fileOffset - file offset to TotalSummaryBlock
   *    isLowToHigh - indicates byte order is low to high if true, else is high to low
   * */
    BBTotalSummaryBlock::BBTotalSummaryBlock(std::ifstream& fis, int64_t fileOffset)
    {
        // save the seekable file handle  and B+ Tree file offset
        this->pfis_ = &fis;
        summaryBlockOffset_ = fileOffset;

        try {
            pfis_->clear();
            pfis_->seekg(fileOffset, std::ios_base::beg);
             if (pfis_->eof())
                std::cerr<<"Hit end of file in seekg in BBTotalSummaryBlock\n";


            pfis_->read( reinterpret_cast<char*>(&basesCovered_) , sizeof(int64_t) );
            basesCovered_=endian::LittleDouble(basesCovered_);

            pfis_->read( reinterpret_cast<char*>(&minVal_) , sizeof(float) );
            minVal_=endian::LittleFloat(minVal_);

            pfis_->read( reinterpret_cast<char*>(&maxVal_) , sizeof(float) );
            maxVal_=endian::LittleFloat(maxVal_);

            pfis_->read( reinterpret_cast<char*>(&sumData_) , sizeof(float) );
            sumData_=endian::LittleFloat(sumData_);

            pfis_->read( reinterpret_cast<char*>(&sumSquares_) , sizeof(float) );
            sumSquares_=endian::LittleFloat(sumSquares_);


        }catch(...) {
            throw std::runtime_error("Error reading Total Summary Block");
            }

        }

    /*
    *   Constructor for filling in TotalSummaryBlock
    * */
     BBTotalSummaryBlock::BBTotalSummaryBlock(int64_t basesCovered, float minVal, float maxVal,
                               float sumData, float sumSquares){

        basesCovered_ = basesCovered;
        minVal_ = minVal;
        maxVal_ = maxVal;
        sumData_ = sumData;
        sumSquares_ = sumSquares;
    }
