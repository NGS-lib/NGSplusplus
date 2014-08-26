#ifndef ZOOMDATABLOCK_H
#define ZOOMDATABLOCK_H


// ***************************************************************************
//   ZoomDataBlock.h (c) 2014
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



#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "RPTreeLeafNodeItem.h"
#include "ZoomDataRecord.h"
#include "endian_helper.h"
#include "decompress_util.h"
#include <stdio.h>
#include <stdlib.h>
//#include

class ZoomDataBlock
{
    public:



       // ZoomDataBlock() {};
        virtual ~ZoomDataBlock() {

        };

        ZoomDataBlock(int32_t zoomLevel, std::ifstream* fis, RPTreeLeafNodeItem* leafHitItem,
                         std::map<uint32_t, std::string>* chromosomeMap, int32_t uncompressBufSize)
                         {

        this->zoomLevel_ = zoomLevel;
        this->leafHitItem_ = leafHitItem;
        this->chromosomeMap_ = chromosomeMap;
        this->fis_=fis;
        fileOffset_ = this->leafHitItem_->getDataOffset();
        dataBlockSize_ = this->leafHitItem_->getDataSize();
        zoomBuffer_.resize(uncompressBufSize);

            //   = new char[(int32_t) dataBlockSize_];
          //  fis->read(&buffer[0],itemBlockSize);
       // std::stringstream       itemStream;
        //Set string stream, does not copy buffer
       // itemStream.rdbuf()->pubsetbuf(&zoomBuffer_[0],dataBlockSize_);

        try {
            fis->clear();
            fis->seekg(fileOffset_);
            if (fis->eof())
                std::cerr<<"Hit end of file in seekg in ZoomDataBlock()\n";


        this->leafHitItem_ = leafHitItem;
        this->chromosomeMap_ = chromosomeMap;

      //  RPTreeLeafNodeItem* itrTest =*leafHitItem_;


        fileOffset_ = this->leafHitItem_->getDataOffset();
        dataBlockSize_ = this->leafHitItem_->getDataSize();
        //char * tempBuffer = new char[(int32_t) leafDataSize_];

        std::istream_iterator<char> streamIt(*fis);
        std::vector<unsigned char> tempBuffer;

        // read Wig data block into a buffer
        for ( int i=0  ; i<dataBlockSize_;  i++)
            {
            tempBuffer.push_back((unsigned char)(*streamIt));
            streamIt++;
            }


/************************************************************
            if (uncompressBufSize > 0)
                read_compressed(*fis_,&zoomBuffer_[0],dataBlockSize_);
            else
                 fis->read(&zoomBuffer_[0],dataBlockSize_);  // use uncompressed read buffer directly

***************************************************************/
        #ifdef VERBOSE
                    std::cerr<<"Was this compressed :"<<(uncompressBufSize>0)<<"\n";
                    std::cerr<<"Buffer resize is "<<dataBlockSize_<<"\n";
                    std::cerr<<"Buffer element count is "<<zoomBuffer_.size()<<"\n";
                    std::cerr<<"Uncompressed Buffer is ";
                    std::copy(tempBuffer.begin(), tempBuffer.end(), std::ostream_iterator<char>(std::cerr));
                    std::cerr<<"\n";
        #endif

        } catch (...) {

           // log.error("Error reading Zoom level " + this.zoomLevel + " data for leaf item ",  ex);
            //String error = String.format("Error reading zoom level %d data for leaf item %d\n", this.zoomLevel);
           // std::string s;
            std::stringstream out;
            out << zoomLevel;
            throw std::runtime_error("Error reading zoom level"+out.str()+" data for leaf item \n");
        }
        // initialize unread data size
        remDataSize_ = dataBlockSize_;

        // use method getZoomData to extract block data
    }

    /*
    *   Method returns all zoom level data within the decompressed block buffer
    *
    *   Parameters:
    *       selectionRegion - chromosome region for selecting zoom level data records
    *       contained - indicates selected data must be contained in selection region
    *           if true, else may int32_tersect selection region
    *
    *   Returns:
    *      zoom data records in the data block
    *
    *   Note: Remaining bytes to data block are used to determine end of reading
    *   since a zoom record count for the data block is not known.
    * */

    std::vector<ZoomDataRecord*> getZoomData(RPChromosomeRegion* selectionRegion,
                                                 bool contained) {

        int32_t chromID, chromStart, chromEnd, validCount;
        float minVal, maxVal, sumData, sumSquares;
        int32_t itemHitValue;
        int32_t recordNumber = 0;

        // allocate the bed feature array list
       // zoomDataList_;

        // check if all leaf items are selection hits
        //TODo VALIDATE stream reading here
        std::stringstream  itemStream;
            //Set string stream, does not copy buffer
        itemStream.rdbuf()->pubsetbuf(&zoomBuffer_[0],dataBlockSize_);

        RPChromosomeRegion* itemRegion = new RPChromosomeRegion(leafHitItem_->getChromosomeBounds());
        int32_t leafHitValue = itemRegion->compareRegions(selectionRegion);

        try {
            //for(int32_t index = 0; mRemDataSize >= ZoomDataRecord.RECORD_SIZE; ++index) {
            for (int32_t index = 0; remDataSize_ > 0; ++index) {
                recordNumber = index + 1;


                    itemStream.read( reinterpret_cast<char*>(&chromID) , sizeof(int32_t) );
                    chromID=endian::LittleLong(chromID);

                    itemStream.read( reinterpret_cast<char*>(&chromStart) , sizeof(int32_t) );
                    chromStart=endian::LittleLong(chromStart);

                    itemStream.read( reinterpret_cast<char*>(&chromEnd) , sizeof(int32_t) );
                    chromEnd=endian::LittleLong(chromEnd);

                    itemStream.read( reinterpret_cast<char*>(&validCount) , sizeof(int32_t) );
                    validCount=endian::LittleLong(validCount);

                    itemStream.read( reinterpret_cast<char*>(&minVal) , sizeof(int32_t) );
                    minVal=endian::LittleFloat(minVal);

                    itemStream.read( reinterpret_cast<char*>(&maxVal) , sizeof(int32_t) );
                    maxVal=endian::LittleFloat(maxVal);

                    itemStream.read( reinterpret_cast<char*>(&sumData) , sizeof(int32_t) );
                    sumData=endian::LittleFloat(sumData);

                    itemStream.read( reinterpret_cast<char*>(&sumSquares) , sizeof(int32_t) );
                    sumSquares=endian::LittleFloat(sumSquares);


                    #ifdef VERBOSE
                    std::cerr<<"value is at ZoomDataLevel"<<chromID<<" "<<chromStart<<" "<<chromEnd<<"\n";
                    #endif


                if (leafHitValue == 0) {     // contained leaf region always a hit
                    std::string chromName = chromosomeMap_->at(chromID);
                    ZoomDataRecord* zoomRecord = new ZoomDataRecord(zoomLevel_, recordNumber, chromName,
                            chromID, chromStart, chromEnd, validCount, minVal, maxVal, sumData, sumSquares);
                    zoomDataList.push_back(zoomRecord);
                } else {      // test for hit
                    itemRegion = new RPChromosomeRegion(chromID, chromStart, chromID, chromEnd);
                    itemHitValue = itemRegion->compareRegions(selectionRegion);

                    // itemHitValue < 2 for int32_tersection; itemHitValue == 0 for is contained
                    if ( ( !contained && std::abs(itemHitValue) < 2) || itemHitValue == 0) {
                        std::string chromName = chromosomeMap_->at(chromID);
                        ZoomDataRecord* zoomRecord = new ZoomDataRecord(zoomLevel_, recordNumber, chromName,
                                chromID, chromStart, chromEnd, validCount, minVal, maxVal, sumData, sumSquares);
                        zoomDataList.push_back(zoomRecord);
                    }
                }

                // compute data block remainder fom size item read
                remDataSize_ -= ZoomDataRecord::RECORD_SIZE;
            }

        } catch (...) {
            // accept this as an end of block condition unless no items were read
            //TODO I think not..do not use throw catch as a flow control mecanism
            if (recordNumber == 1)
                throw std::runtime_error("Read error for zoom level ");
        }

        return zoomDataList;
    }

    protected:
    private:


    // Bed data block access variables   - for reading in bed records from a file
    int64_t fileOffset_;       // data block file offset
    int64_t dataBlockSize_;    // byte size for data block specified in the R+ leaf
    std::ifstream* fis_;
    // defines the zoom level source chromosomes
    int32_t zoomLevel_;         // zoom level for the R+ chromosome data location tree
    std::map<uint32_t, std::string>* chromosomeMap_;  // map of chromosome ID's and corresponding names
    RPTreeLeafNodeItem* leafHitItem_;   //R+ leaf item with chromosome region and file data location

    // Provides uncompressed byte stream data reader
    std::vector<char>  zoomBuffer_;  // buffer containing leaf block data uncompressed
    int32_t remDataSize_;   // number of unread decompressed data bytes

    // Bed data extraction members
    std::vector<ZoomDataRecord*> zoomDataList; // array of zoom level data

};




#endif // ZOOMDATABLOCK_H
