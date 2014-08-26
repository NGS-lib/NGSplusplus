#ifndef BBZOOMLEVELFORMAT_H
#define BBZOOMLEVELFORMAT_H

// ***************************************************************************
//    BBZoomLevelFormat.h (c) 2014
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

#include <string>
#include <iostream>
#include "stdint.h"
#include <stdexcept>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
class BBZoomLevelFormat
{
    public:
        BBZoomLevelFormat();
        virtual ~BBZoomLevelFormat();


         BBZoomLevelFormat(uint32_t zoomLevel, std::ifstream& fis, uint64_t fileOffset,
                           uint64_t dataSize, uint32_t uncompressBufSize);

        uint32_t getZoomLevel_() { return zoomLevel_; }
        //void SetZoomLevel_(uint32_t val) { zoomLevel_ = val; }
        uint64_t getZoomFormatOffset() { return zoomFormatOffset_; }
        //void SetZoomFormatOffset(uint64_t val) { zoomFormatOffset_ = val; }
        uint64_t getZoomDataOffset() { return zoomDataOffset_; }
        //void SetZoomDataOffset(uint64_t val) { zoomDataOffset_ = val; }
        uint64_t getZoomIndexOffset() { return zoomIndexOffset_; }
        //void SetZoomIndexOffset(uint64_t val) { zoomIndexOffset_ = val; }
        uint32_t getZoomRecordCount() { return zoomRecordCount_; }
        //void SetZoomRecordCount(uint32_t val) { zoomRecordCount_ = val; }
        uint64_t GetZoomDataSize() { return zoomDataSize_; }
        //void SetZoomDataSize(uint64_t val) { zoomDataSize_ = val; }

        const static uint32_t ZOOM_FORMAT_HEADER_SIZE = 4;
        const  static uint32_t MAX_ZOOM_DATA_RECORDS = 100000000;



    protected:
    private:


        uint32_t zoomLevel_;
        std::ifstream* pfis_;
        uint64_t zoomFormatOffset_;
        uint64_t zoomDataOffset_;
        uint64_t zoomIndexOffset_;
        uint32_t zoomRecordCount_;
        uint64_t zoomDataSize_;
};


#endif // BBZOOMLEVELFORMAT_H
