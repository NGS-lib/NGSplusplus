
// ***************************************************************************
//   BBZoomLevels.cpp (c) 2014
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


#include "BBZoomLevels.h"
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "stdint.h"

BBZoomLevels::BBZoomLevels()
{
    //ctor
}

BBZoomLevels::~BBZoomLevels()
{
    //dtor
}


  /*
   *  constructor   - reads zoom level headers and data format from file I/O stream
   *      for all zoom levels and acts as a container for zoom data records.
   *
   *  Parameters:
   *      fis - file input stream handle
   *      fileOffset - file byte location for zoom level headers
   *      zoomLevels - count of zoom levels from BBFile Table C
   *      isLowToHigh - bool flag indicates if values are arranged low to high bytes.
   *      uncompressBufSize - byte size of the buffer to use for decompression
   * */

     BBZoomLevels::BBZoomLevels(std::ifstream& fis, uint64_t fileOffset, uint32_t zoomLevels
                        , uint32_t uncompressBufSize){
        uint32_t zoomLevel;
        uint32_t zoomHeadersRead;
        uint64_t zoomDataOffset;
        uint64_t zoomIndexOffset;

        // save the seekable file handle and zoom zoomLevel headers file offset
        //this.fis = fis;à

        zoomHeadersOffset_ = fileOffset;
        zoomLevelsCount_ = zoomLevels;

        // Note: a bad zoom header will result in a 0 count returned
        zoomHeadersRead =  readZoomHeaders(fis, zoomHeadersOffset_, zoomLevels);

        if(zoomHeadersRead > 0){

            // create zoom level data format containers
//            zoomLevelFormatList = new std::vector<BBZoomLevelFormat>();

            // for each zoom zoomLevel, get associated zoom data format


            for(uint32_t index = 0; index < zoomHeadersRead; ++index) {

                zoomLevel = index + 1;

                // Zoom dataOffset (from Table D) is file location for zoomCount (Table O)
                // Note: This dataOffset is zoomFormatLocation in BBZoomLevelFormat.
                zoomDataOffset = zoomLevelHeaders_.at(index)->getDataOffset();
         //       uint64_t indexoffset = zoomLevelHeaders_.at(index)->getIndexOffset() ;
                // R+ zoom index offset (Table D) marks end of zoom data in the
                // zoom level format (Table O)
                uint64_t dataSize = zoomLevelHeaders_.at(index)->getIndexOffset() - zoomDataOffset
                        - BBZoomLevelFormat::ZOOM_FORMAT_HEADER_SIZE;

                // get zoom zoomLevel data records  - zoomDataOffset references zoomCount in Table O
                // Note: zoom zoomLevel data records read their own data
                BBZoomLevelFormat* zoomLevelData = new BBZoomLevelFormat(zoomLevel, fis, zoomDataOffset,
                        dataSize, uncompressBufSize);

                zoomLevelFormatList_.push_back(zoomLevelData);
            }

            // for each zoom zoomLevel, get associated R+ tree
            for(uint32_t index = 0; index < zoomHeadersRead; ++index) {

                // Zoom indexOffset (from Table D) is file location
                // for Table O zoomIndex for R+ tree zoom data
                zoomIndexOffset = zoomLevelHeaders_.at(index)->getIndexOffset();

                // get Zoom Data R+ Tree (Tables K, L, M, N): exists for zoom levels
                RPTree* zoomRPTree = new RPTree(fis, zoomIndexOffset, uncompressBufSize, true);

                //if(zoomRPTree.getNodeCount() > 0)
                    zoomLevelRPTree_.push_back(zoomRPTree);
            }
        }
    }

    /*
    *   Method returns the BBFile's first zoom header file offset.
    *
    *   Note zoom headers immediately follow the BBFile header (Table C)
    *
    *   Returns:
    *       first zoom header file offset
    * */
    uint64_t BBZoomLevels::getZoomHeadersOffset() {
        return zoomHeadersOffset_;
    }

    /*
    *   Method returns the number of zoom level headers found.
    *
    *   Note Should match zoomLevels in the BBFile header (Table C)
    *
    *   Returns:
    *      number of zoom level headers found
    * */
    uint32_t BBZoomLevels::getZoomHeaderCount() {
        return zoomLevelHeaders_.size();
    }

    /*
    *   Method returns the zoom level headers.
    *
    *   Returns:
    *      zoom level headers
    * */
     std::vector<BBZoomLevelHeader*> BBZoomLevels::getZoomLevelHeaders() {
        return zoomLevelHeaders_;
    }

    /*
    *   Method returns the zoom level header for specified level.
    *
    *   Parameters:
    *       level - zoom level; level starts at 1
    *
    *   Returns:
    *      Zoom level header for specified level; or null for bad zoom level.
    * */
    BBZoomLevelHeader* BBZoomLevels::getZoomLevelHeader(uint32_t level) {
        if(level < 1 || level > zoomHeadersOffset_)
            return NULL;

        return zoomLevelHeaders_.at(level - 1);
    }

    /*
    *   Method returns the zoom level formats for zoom data.
    *
    *   Returns:
    *      zoom level formats for zoom data
    * */
    std::vector<BBZoomLevelFormat*> BBZoomLevels::getZoomLevelFormats(){
        return zoomLevelFormatList_;
    }

    /*
    *   Method returns the R+ index tree for the specified zoom level.
    *
    *   Parameters:
    *       level - zoom level; level starts at 1
    *
    *   Returns:
    *      R+ index tree for the specified zoom level; or null for bad zoom level
    * */
    RPTree* BBZoomLevels::getZoomLevelRPTree(uint32_t level) {
        if(level < 1 || level > zoomHeadersOffset_)
            return NULL;

        return zoomLevelRPTree_.at(level - 1);
    }

    /*
    * Reads in all the Zoom Headers.
    *
    *   Parameters:
    *       fileOffset - File byte location for first zoom level header
    *       zoomLevels - count of zoom levels to read in
    *       isLowToHigh - indicate byte order is lwo to high, else is high to low
    *
    *   Returns:
    *       Count of zoom levels headers read, or 0 for failure to find the
    *       header information.
    * */
    //private
    uint32_t BBZoomLevels::readZoomHeaders(std::ifstream& fis, uint64_t fileOffset, uint32_t zoomLevels) {
        uint32_t level = 0;
        BBZoomLevelHeader* zoomLevelHeader;

        if(zoomLevels < 1)
            return 0;

        // get zoom header information for each zoom levelsRead
        for(uint32_t index = 0; index < zoomLevels; ++index)  {
            level = index + 1;

            // read zoom level header - read error is returned as Runtime Exception
            zoomLevelHeader = new BBZoomLevelHeader(fis, fileOffset, level);

            zoomLevelHeaders_.push_back(zoomLevelHeader);

            fileOffset += BBZoomLevelHeader::ZOOM_LEVEL_HEADER_SIZE;
        }

        return level;
    }
