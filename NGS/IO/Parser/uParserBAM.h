// ***************************************************************************
// uParserBAM.h (c) 2013
// Alexei Nordell-Markovits : Sherbrooke University
// Charles Joly Beauparlant : Laval University
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



#ifndef UPARSEBAM_H_INCLUDED
#define UPARSEBAM_H_INCLUDED
#include "uParserBase.h"
#include <iostream>
#include "api/BamReader.h"
namespace NGS
{

class uParserBAM: public uParserBase
{

public :
    uParserBAM();
    ~uParserBAM();
    void init(const std::string& filename, bool header = false);
    void init(std::istream* stream, bool header = false);

    uToken getNextEntry();
	static uParserBase * Create() { return new uParserBAM(); }
	bool eof();
private:
 //   static DerivedParserRegister<uParserBAM> reg;
    BamTools::BamReader m_BamReader;
    bool m_IsBuffer=false;
    BamTools::BamAlignment m_BufferAlignement;
    void _parseHeader();
    /**< String for dynamic parsing of Sam */
};

}


#endif // UPARSEBAM_H_INCLUDED
