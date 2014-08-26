#ifndef BBZOOMLEVELHEADER_H
#define BBZOOMLEVELHEADER_H

// ***************************************************************************
//   BBZoomLevelHeader.h (c) 2014
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

#include <vector>
#include "stdint.h"
#include <string>
class BBZoomLevelHeader
{
    typedef char byte;
    private:
        uint64_t zoomLevelheaderOffset_;
        uint32_t zoomLevel_;
        uint32_t reductionLevel_;
        uint32_t reserved_;
        uint64_t dataOffset_;
        uint64_t indexOffset_;

        std::ifstream* pfis_;          // BBFile input stream handle

        void readZoomLevelHeader(uint64_t fileOffset, uint32_t zoomLevel);

    public:
        BBZoomLevelHeader();
        ~BBZoomLevelHeader();

        static const uint32_t ZOOM_LEVEL_HEADER_SIZE = 24;
        uint64_t getZoomLevelheaderOffset() { return zoomLevelheaderOffset_; }
        void setZoomLevelheaderOffset(uint64_t val) { zoomLevelheaderOffset_ = val; }
        uint32_t getZoomLevel() { return zoomLevel_; }
        void setZoomLevel(uint32_t val) { zoomLevel_ = val; }
        uint32_t getReductionLevel() { return reductionLevel_; }
        void setReductionLevel(uint32_t val) { reductionLevel_ = val; }
        uint32_t getReserved() { return reserved_; }
        void setReserved(uint32_t val) { reserved_ = val; }
        uint64_t getDataOffset() { return dataOffset_; }
        void setDataOffset(uint64_t val) { dataOffset_ = val; }
        uint64_t getIndexOffset() { return indexOffset_; }
        void setIndexOffset(uint64_t val) { indexOffset_ = val; }

     BBZoomLevelHeader(std::ifstream& fis, uint64_t fileOffset, uint32_t zoomLevel);
     BBZoomLevelHeader(uint32_t zoomLevel, uint32_t reductionLevel, uint32_t reserved,
                             uint64_t dataOffset, uint64_t indexOffset);


};

#endif // BBZOOMLEVELHEADER_H
