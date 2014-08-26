#ifndef BBTOTALSUMMARYBLOCK_H
#define BBTOTALSUMMARYBLOCK_H

// ***************************************************************************
//   BBTotalSummaryBlock.h (c) 2014
//   Copyright @ Alexei Nordell-Markovits : Sherbrooke University
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

#include <string>
#include <iostream>
#include "stdint.h"
#include <stdexcept>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>

class BBTotalSummaryBlock
{


    private:
        std::ifstream* pfis_;
        uint64_t summaryBlockOffset_;
        uint64_t basesCovered_;
        float minVal_;
        float maxVal_;
        float sumData_;
        float sumSquares_;

    public:

        BBTotalSummaryBlock(std::ifstream& fis, int64_t fileOffset);
        BBTotalSummaryBlock(int64_t basesCovered, float minVal, float maxVal,
                                   float sumData, float sumSquares);

        uint64_t GetSummaryBlockOffset_() { return summaryBlockOffset_; }
        void SetSummaryBlockOffset_(int64_t val) { summaryBlockOffset_ = val; }
        uint64_t GetBasesCovered_() { return basesCovered_; }
        void SetBasesCovered_(int64_t val) { basesCovered_ = val; }
        float GetMinVal_() { return minVal_; }
        void SetMinVal_(float val) { minVal_ = val; }
        float GetMaxVal_() { return maxVal_; }
        void SetMaxVal_(float val) { maxVal_ = val; }
        float GetSumData_() { return sumData_; }
        void SetSumData_(float val) { sumData_ = val; }
        float GetSumSquares_() { return sumSquares_; }
        void SetSumSquares_(float val) { sumSquares_ = val; }

        static const uint32_t TOTAL_SUMMARY_BLOCK_SIZE = 40;

    protected:





};

#endif // BBTOTALSUMMARYBLOCK_H
