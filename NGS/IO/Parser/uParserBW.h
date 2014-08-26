#ifndef UPARSERBW_H_INCLUDED
#define UPARSERBW_H_INCLUDED

// ***************************************************************************
// uParserBW.h (c) 2014
// Alexei Nordell-Markovits : Sherbrooke University
//
//       This file is part of the NGS++ library.
//
//    The NGS++ library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program (lgpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// ***************************************************************************
#include "uParserBase.h"
#include "../../utility/utility.h"
#include <iostream>
#include "../../bbfile/include/BBFileReader.h"
//#include "uParserFactory.h"
namespace NGS
{

class uParserBW: public uParserBase
{

public :
    uParserBW();
    ~uParserBW();
    virtual void init(const std::string& filename, bool header = false);
    virtual void init(std::istream* stream, bool header = false);

    uToken getNextEntry();
	static uParserBase * Create() { return new uParserBW(); }
	bool eof();
private:
 //   static DerivedParserRegister<uParserBAM> reg;
    BBFileReader m_BBReader;
    BigWigIterator m_BBIterator;
    std::ifstream m_fstream;
    bool m_IsBuffer=false;
};

}
#endif // UPARSERBW_H_INCLUDED
