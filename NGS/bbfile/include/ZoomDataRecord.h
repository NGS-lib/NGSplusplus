#ifndef ZOOMDATARECORD_H
#define ZOOMDATARECORD_H


// ***************************************************************************
//   ZoomDataRecord.h (c) 2014
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
#include "stdint.h"

class ZoomDataRecord
{
    public:
        ZoomDataRecord() {}
        virtual ~ZoomDataRecord() {}

        const static int32_t RECORD_SIZE = 32;

   /*
    *   Constructor for filling in zoom data record class.
    *
    *   Parameters:
    *       zoomLevel - level of zoom
    *       recordNumber - record sequence number of multiple zoom level records
    *       chromName - chromosome/contig name
    *       chromId - mChromosome ID
    *       chromstart - starting base for zoom data region
    *       chromEnd - ending base for zoom data region
     *      validCount - number of bases in the region for which there is data
     *      minVal - minimum value in region
     *      maxVal - maximum value in region
     *      sumData - sum of all region data
     *      sumSquares - sum of the squares of all region data
     *
    * */
     ZoomDataRecord(int32_t zoomLevel, int32_t recordNumber, std::string chromName, int32_t chromId, int32_t chromStart, int32_t chromEnd,
            int32_t validCount, float minVal, float maxVal, float sumData, float sumSquares ){

        this->zoomLevel_ = zoomLevel;
        this->recordNumber_ = recordNumber;
        this->chromName_ = chromName;
        this->chromId_ = chromId;
        this->chromStart_ = chromStart;
        this->chromEnd_ = chromEnd;
        this->basesCovered_ = validCount;
        this->minVal_ = minVal;
        this->maxVal_ = maxVal;
        this->sumData_ = sumData;
        this->sumSquares_ = sumSquares;
    }

    int32_t getZoomLevel() {
        return zoomLevel_;
    }

    int32_t getRecordNumber() {
        return recordNumber_;
    }

     std::string getChromName() {
        return chromName_;
    }

     int32_t getChromId() {
        return chromId_;
    }

    int32_t getChromStart() {
        return chromStart_;
    }

    int32_t getChromEnd() {
        return chromEnd_;
    }

    int32_t getBasesCovered() {
        return basesCovered_;
    }

    float getMinVal() {
        return minVal_;
    }

    float getMaxVal() {
        return maxVal_;
    }

    float getSumData() {
        return sumData_;
    }

    float getMeanVal() {
        return basesCovered_ == 0 ? 0 : sumData_ / basesCovered_;
    }

    float getSumSquares() {
        return sumSquares_;
    }

    protected:
    private:

    int32_t zoomLevel_;         // zoom level associated with data
    int32_t recordNumber_;      // record number

    // chromosome region statistics (useful for calculating mean and standard deviation)
    std::string chromName_;      // chromosome/contig name
    int32_t chromId_;           // Numerical ID for mChromosome/contig
    int32_t chromStart_;        // starting base position  (from 0)
    int32_t chromEnd_;          // ending base position
    int32_t basesCovered_;        // number of bases with data
    float minVal_;          // minimum value for file data
    float maxVal_;          // maximum value for file data
    float sumData_;         // sum of all squares of file data values
    float sumSquares_;      // sum of squares of file data values
};

#endif // ZOOMDATARECORD_H
