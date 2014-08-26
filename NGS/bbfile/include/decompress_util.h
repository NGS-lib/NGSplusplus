#ifndef DECOMPRESS_UTIL_H_INCLUDED
#define DECOMPRESS_UTIL_H_INCLUDED


// ***************************************************************************
//   decompress_util.h (c) 2014
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



#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include "endian_helper.h"
#include <vector>

using namespace boost::iostreams;


static void return_compressed(std::ifstream&stream,char* buffer, int32_t read_size)
{

 // std::ifstream datfile(filename, std::ios::binary);
  boost::iostreams::filtering_istreambuf zdat;
  zdat.push(boost::iostreams::zlib_decompressor());
  zdat.push(stream);
  boost::iostreams::read(zdat, buffer, read_size);
}

static void read_compressed(std::ifstream& stream, char* buffer, int32_t read_size){

    try {
    if ((stream.bad()))
        throw std::runtime_error("Invalid stream object in read_compressed");

  //  boost::iostreams::filtering_stream<boost::iostreams::input> decompressor;
    stream.clear();
    filtering_stream<input> decompressor;
    decompressor.push(zlib_decompressor());
    decompressor.push(stream);

    copy(decompressor, buffer);

     if (decompressor.good()==false){
        std::cerr<<"Stream is bad after read\n";
    }

    }catch(...){throw;}

}

static void gunzip_buffer(std::vector<char>& buffer){

/*
boost::iostreams::filtering_istreambuf in;
in.push(boost::iostreams::zlib_decompressor());
in.push(boost::interprocess::basic_vectorstream<std::vector<char>>(buffer));

std::vector<char> decomp;
boost::iostreams::copy(in, boost::interprocess::basic_vectorstream < std::vector < char >> (decomp));


*/
    std::stringstream  itemStream;
//Set string stream, does not copy buffer
    itemStream.rdbuf()->pubsetbuf(&buffer[0],buffer.size());

    filtering_stream<input> decompressor;
    decompressor.push(zlib_decompressor());
    decompressor.push(itemStream);
    std::stringstream ss;
    copy(decompressor, std::cout);

}

#endif // DECOMPRESS_UTIL_H_INCLUDED
