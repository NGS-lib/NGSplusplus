// ***************************************************************************
// uParserBed.h (c) 2013
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


#ifndef UPARSERBED_H_INCLUDED
#define UPARSERBED_H_INCLUDED

#include "uParserBase.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>
namespace NGS
{

class uParserBed : public uParserBase
{
public :
    uParserBed();
    ~uParserBed();

    virtual void init(const std::string& filename, bool header = false);
    virtual void init(std::istream* stream, bool header = false);
    static uParserBase * Create() { return new uParserBed(); }
    uToken getNextEntry();

protected:
    char m_delimiter = '\t';

private:
    void _parseHeader();
    int _countColumns(char* line) const;
    void _validateColumnNumber(int numberOfColumn) const;
    void _convertLineToTokenInfosBed(char* line, std::stringstream& token_infos);
    void _pushBackLine(char* line);
    std::string _getNextEntry(char* line);

    int m_numberOfColumn = 0;
    bool m_headerParsed = false;
  //  static DerivedParserRegister<uParserBed> reg;
};

} // End of namespace NGS
#endif // UPARSERBED_H_INCLUDED
