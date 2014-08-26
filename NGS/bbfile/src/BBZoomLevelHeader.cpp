
// ***************************************************************************
//   BBZoomLevelHeader.cpp (c) 2014
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


#include "BBZoomLevelHeader.h"
#include <iostream>
#include <fstream>
#include "stdint.h"
#include <stdexcept>
#include "endian_helper.h"

BBZoomLevelHeader::BBZoomLevelHeader()
{
    //ctor
}

BBZoomLevelHeader::~BBZoomLevelHeader()
{
    //dtor
}

 /*
    *   Constructor reads zoom level header
    *
    *   Parameters:
    *       fis - File input stream handle
    *       fileOffset - file byte position for zoom header
    *       zoomLevel - level of zoom
    * */
   BBZoomLevelHeader::BBZoomLevelHeader(std::ifstream& fis, uint64_t fileOffset, uint32_t zoomLevel){

        this->pfis_ = &fis;
        zoomLevelheaderOffset_ = fileOffset;
        this->zoomLevel_ = zoomLevel;

        readZoomLevelHeader(zoomLevelheaderOffset_, this->zoomLevel_);
    }

    /*
    *   Constructor loads zoom level header according to parameter specification.
    *
    *   Parameters: (as defined above)
    * */
    BBZoomLevelHeader::BBZoomLevelHeader(uint32_t zoomLevel, uint32_t reductionLevel, uint32_t reserved,
                             uint64_t dataOffset, uint64_t indexOffset){
        zoomLevel_ = zoomLevel;
        reductionLevel_ = reductionLevel;
        reserved_ = reserved;
        dataOffset_ = dataOffset;
        indexOffset_ = indexOffset;
    }

  /*
    *   Reads zoom level header information into class data members.
    *
    *   Parameters:
    *       fileOffset - Byte position in fle for zoom header
    *       zoomLevel - level of zoom
    * */
    void BBZoomLevelHeader::readZoomLevelHeader(uint64_t fileOffset, uint32_t zoomLevel) {
    typedef char byte;
      // LittleEndianInputStream lbdis = null;
      // DataInputStream bdis = null;
            try {

            // Read zoom header into a buffer
            pfis_->clear();
            pfis_->seekg(fileOffset, std::ios_base::beg);
            if (pfis_->eof())
                std::cerr<<"Hit end of file in seekg in ReadZoomLevelHeader\n";
            // Get zoom level information
            pfis_->read( reinterpret_cast<char*>(&reductionLevel_) , sizeof(uint32_t) );
            reductionLevel_=endian::LittleLong(reductionLevel_);

            pfis_->read( reinterpret_cast<char*>(&reserved_) , sizeof(uint32_t) );
            reserved_=endian::LittleLong(reserved_);

            pfis_->read( reinterpret_cast<char*>(&dataOffset_) , sizeof(uint64_t) );
            dataOffset_=endian::LittleDouble(dataOffset_);

            pfis_->read( reinterpret_cast<char*>(&indexOffset_) , sizeof(uint64_t) );
            indexOffset_=endian::LittleDouble(indexOffset_);


        }catch(...) {
            throw std::runtime_error("Error reading zoom header");
        }
    }
