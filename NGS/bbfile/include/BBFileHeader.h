#ifndef BWFILEHEADER_H_INCLUDED
#define BWFILEHEADER_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>

#include "stdint.h"
#include "endian_helper.h"
#include <stdexcept>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>


// ***************************************************************************
//   BBFileHeader.h (c) 2014
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





/**
 *  : Alexei Nordell Markovits
 * Refactor of IGV BigWig Java reader
 */

/**
 *  Container class defines the header information for BigBed and BigWig files
 */



class BBFileHeader {
    typedef unsigned char byte;
    private:
    //What is this
  //  static Logger log = Logger.getLogger(BBFileHeader.class);

    public:
    // defines bigBed/bigwig Header Format types
    static const uint32_t BBFILE_HEADER_SIZE = 64;

    static const uint32_t BIGWIG_MAGIC_LTH = 0x888FFC26; // BigWig Magic Low to High
    static const uint32_t BIGWIG_MAGIC_HTL = 0x26FC8F66; // BigWig Magic High to Low

    static const uint32_t BIGBED_MAGIC_LTH = 0x8789F2EB; // BigBed Magic Low to High
    static const uint32_t BIGBED_MAGIC_HTL = 0xEBF28987; // BigBed Magic High to Low

    private :
    // defines the bigBed/bigWig source file access
     std::string path;               // bigBed file/pathname
     std::ifstream* fis;      // BBFile I/O stream handle
      //std::ifstream ifs("test.txt", ios::in | ios::binary);
     uint64_t fileHeaderOffset;     // file offset for file header

     bool isHeaderOk_;        // File header read correctly?
    // bool isLowToHigh_;       // flag indicates values represented low to high bytes
     bool isBigBed_;          // flag indicates file is BigBed format
     bool isBigWig_;          // flag indicates file is BigWig format;

    // BBFile Header items - Table C:
    // mMagic number (4 bytes) indicates file type and byte order :
    // 0x888FFC26 for bigWig, little endian if swapped
    // 0x8789F2EB for bigBed, little endian if swapped
    uint32_t magic_;                // 4 byte mMagic Number
    uint16_t version_;            // 2 byte version ID; currently 3
    uint16_t nZoomLevels_;         // 2 byte count of zoom sumary levels
    uint64_t chromTreeOffset_;     // 8 byte offset to mChromosome B+ Tree index
    uint64_t fullDataOffset_;      // 8 byte offset to unzoomed data dataCount
    uint64_t fullIndexOffset_;     // 8 byte offset to R+ Tree index of items
    uint16_t fieldCount_;         // 2 byte number of fields in bed. (0 for bigWig)
    uint16_t definedFieldCount_;  // 2 byte number of fields that are bed fields
    uint64_t autoSqlOffset_;       // 8 byte offset to 0 terminated string with .as spec
    uint64_t totalSummaryOffset_;  // 8 byte offset to file summary data block
    uint32_t uncompressBuffSize_;  // 4 byte maximum size for decompressed buffer
    uint64_t reserved_;            // 8 bytes reserved for future expansion. Currently 0

    public:
    // constructor reads BBFile header from an input stream
    BBFileHeader(std::string path, int64_t fileOffset) {
        try {
            // save the path and seekable file handle
            this->path = path;

            this->fis= new std::ifstream(path.c_str(), std::ifstream::in | std::ifstream::binary);
         //   this->fis->open();
            fileHeaderOffset = fileOffset;
            // read in BBFile header
            this->isHeaderOk_ = readBBFileHeader(fileHeaderOffset);
        }
        catch(...){
            throw;
        }
    }

    BBFileHeader(std::string path, std::ifstream& fis, int64_t fileOffset) {

        // save the path and seekable file handle
        this->path = path;
        this->fis = &fis;
        fileHeaderOffset = fileOffset;

        // read in BBFile header
        isHeaderOk_ = readBBFileHeader(fileHeaderOffset);

    }


    /*
    *   Constructor loads BBFile header class from parameter specifications.
    *
    *   Parameters: (as defined above)
    * */
    BBFileHeader(
            uint32_t magic,
            uint16_t version,
            uint16_t zoomLevels,
            uint64_t chromTreeOffset,
            uint64_t fullDataOffset,
            uint64_t fullIndexOffset,
            uint16_t fieldCount,
            uint16_t definedFieldCount,
            uint64_t autoSqlOffset,
            uint64_t totalSummaryOffset,
            uint32_t uncompressBuffSize,
            uint64_t reserved) {

        this->magic_ = magic;
        // Note: may want to validate the rest of the fields as well
        if (isBigWig() || isBigBed())
            this->isHeaderOk_ = true;

        this->version_ = version;
        this->nZoomLevels_ = zoomLevels;
        this->chromTreeOffset_ = chromTreeOffset;
        this->fullDataOffset_ = fullDataOffset;
        this->fullIndexOffset_ = fullIndexOffset;
        this->fieldCount_ = fieldCount;
        this->definedFieldCount_ = definedFieldCount;
        this->autoSqlOffset_ = autoSqlOffset;
        this->totalSummaryOffset_ = totalSummaryOffset;
        this->uncompressBuffSize_ = uncompressBuffSize;
        this->uncompressBuffSize_ = uncompressBuffSize;
        this->reserved_ = reserved;
    }

     std::string getPath() {
        return path;
    }


    bool isHeaderOK() {
        return isHeaderOk_;
    }

    bool isBigBed() {
        return isBigBed_;
    }

    bool isBigWig() {
        return isBigWig_;
    }

    uint32_t getFileHeaderSize() {
        return BBFILE_HEADER_SIZE;
    }

    // ************* return header items ****************
    uint32_t getMagic() {
        return magic_;
    }

    uint16_t getVersion() {
        return version_;
    }

    uint16_t getZoomLevels() {
        return nZoomLevels_;
    }

    uint64_t getChromosomeTreeOffset() {
        return chromTreeOffset_;
    }

    uint64_t getFullDataOffset() {
        return fullDataOffset_;
    }

    uint64_t getFullIndexOffset() {
        return fullIndexOffset_;
    }

     uint16_t getFieldCount() {
         return fieldCount_;
     }

     uint16_t getDefinedFieldCount() {
         return definedFieldCount_;
     }

    uint64_t getAutoSqlOffset() {
        return autoSqlOffset_;
    }

    uint64_t getTotalSummaryOffset() {
        return totalSummaryOffset_;
    }

    uint32_t getUncompressBuffSize() {
        return uncompressBuffSize_;
    }

    /*
     *  Reads in BBFile header information.
     *
     *  Returns:
     *      Success status flag is true for successfully read header,
     *      or is false for a read error.
    **/
    private :
    bool readBBFileHeader(int64_t fileOffset) {

       // Skipe offset if any
        if (fileOffset>0){
            fis->seekg(fileOffset, std::ios_base::beg);
            if (fis->eof())
                std::cerr<<"Hit end of file in seekg in readBBFileHeader\n";
        }
       // LittleEndianInputStream lbdis = NULL;
      //  DataInputStream bdis = NULL;

       // byte[] buffer = new byte[BBFILE_HEADER_SIZE];

        //Validate file
        try {
            // convert everything to LowToHigh
            // first assume byte order is low to high
            uint32_t header;
            fis->read( reinterpret_cast<char*>(&header) , sizeof(uint32_t) );
            magic_=endian::LittleLong(header);
            // check for a valid bigBed or bigWig file
            if (magic_ == BIGWIG_MAGIC_LTH)
                isBigWig_ = true;
            else if (magic_ == BIGBED_MAGIC_LTH)
                isBigBed_ = true;
            else
                return false;   // can't identify BBFile type

            fis->read( reinterpret_cast<char*>(&version_) , sizeof(uint16_t) );
            version_ =endian::LittleShort(version_);

            fis->read( reinterpret_cast<char*>(&nZoomLevels_) , sizeof(uint16_t) );
            nZoomLevels_ =endian::LittleShort(nZoomLevels_);

            fis->read( reinterpret_cast<char*>(&chromTreeOffset_) , sizeof(uint64_t) );
            chromTreeOffset_ = endian::LittleDouble(chromTreeOffset_);

            fis->read( reinterpret_cast<char*>(&fullDataOffset_) , sizeof(uint64_t) );
            fullDataOffset_ = endian::LittleDouble(fullDataOffset_);

            fis->read( reinterpret_cast<char*>(&fullIndexOffset_) , sizeof(uint64_t) );
            fullIndexOffset_ = endian::LittleDouble(fullIndexOffset_);

            fis->read( reinterpret_cast<char*>(&fieldCount_) , sizeof(uint16_t) );
            fieldCount_ = endian::LittleShort(fieldCount_);

            fis->read( reinterpret_cast<char*>(&definedFieldCount_) , sizeof(uint16_t) );
            definedFieldCount_ = endian::LittleShort(definedFieldCount_);

            fis->read( reinterpret_cast<char*>(&autoSqlOffset_) , sizeof(uint64_t) );
            autoSqlOffset_ = endian::LittleDouble(autoSqlOffset_);

            fis->read( reinterpret_cast<char*>(&totalSummaryOffset_) , sizeof(uint64_t) );
            totalSummaryOffset_ = endian::LittleDouble(totalSummaryOffset_);

            fis->read( reinterpret_cast<char*>(&uncompressBuffSize_) , sizeof(uint32_t) );
            uncompressBuffSize_ = endian::LittleLong(uncompressBuffSize_);

            fis->read( reinterpret_cast<char*>(&reserved_) , sizeof(uint64_t) );
            reserved_ = endian::LittleDouble(reserved_);

        } catch (...) {
            throw std::runtime_error("Error reading file header for " + path);
        }

        // file header was read properly
        return true;
    }

};  // mEndBase of class BBFileHeader

#endif // BWFILEHEADER_H_INCLUDED
