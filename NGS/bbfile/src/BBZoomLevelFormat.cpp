
// ***************************************************************************
//   BBZoomLevelFormat.cpp (c) 2014
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


#include "BBZoomLevelFormat.h"
#include "endian_helper.h"
#include <iostream>
#include <fstream>
#include "stdint.h"
#include <stdexcept>

BBZoomLevelFormat::BBZoomLevelFormat()
{
    //ctor
}

BBZoomLevelFormat::~BBZoomLevelFormat()
{
    //dtor
}

 /*
    *   constructor   - reads zoom level format (but not the data) for a zoom level.
    *
    *   Parameters:
    *       zoomLevel - zoom level; from 1 to zoomLevels in Table C File Header
    *       fis - input stream reader handle
    *       fileOffset - file location for the zoom format table
    * `     datasize - byte size of (compressed/uncompressed) zoom data block
    *       isLowToHigh - boolean flag indicates if buffer data is arranged low to high byte
    *       uncompressBufSize - buffer size for decompressed data; or 0 for uncompressed
    *
    *   Note: Zoom level data Table O is arranged as:
    *       zoomCount - 4 bytes
    *       zoomData - dataSize bytes
    *       R+ zoom index starts at fileOffset + dataSize
    * */
    BBZoomLevelFormat::BBZoomLevelFormat(uint32_t zoomLevel, std::ifstream& fis, uint64_t fileOffset,
                           uint64_t dataSize, uint32_t uncompressBufSize) {

        // store file access info
        this->zoomLevel_ = zoomLevel;
        this->pfis_ = &fis;
        zoomFormatOffset_ = fileOffset;
        zoomDataSize_ = dataSize;

        // Note: a bad zoom data header will result in a 0 count returned
        // or an IOException

        // size of buffer is ZOOM_FORMAT_HEADER_SIZE to get the record count
        try {

            // Read zoom level data format into a buffer
            pfis_->clear();
            pfis_->seekg(zoomFormatOffset_, std::ios_base::beg);
            if (pfis_->eof())
                std::cerr<<"Hit end of file in seekg in BBZoomLevelFormat\n";
            // dec  ode header - or fail
            pfis_->read( reinterpret_cast<char*>(&zoomRecordCount_) , sizeof(uint32_t) );
            zoomRecordCount_=endian::LittleLong(zoomRecordCount_);


        } catch (...) {
            throw std::runtime_error("Error reading zoom level data records (Table O)");
            }

        // integrity check - should be > 0 or less than a max like 100M records?
            // Note: if trouble reading zoom data records, readAllZoomLevelRecords returns 0
            if(zoomRecordCount_ < 0 || zoomRecordCount_ > MAX_ZOOM_DATA_RECORDS)
                return;  // terminate if bad zoom level data encountered

            // Position file offset past the current zoom level header to pick up
            // the zoom data records which immediately follow.
            zoomDataOffset_ = zoomFormatOffset_ + ZOOM_FORMAT_HEADER_SIZE;

            // calculate the position of the R+ zoom index tree
            zoomIndexOffset_ = zoomDataOffset_ + zoomDataSize_;
    }

