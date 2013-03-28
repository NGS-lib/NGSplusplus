// ***************************************************************************
// uParserBedGraph.h (c) 2013
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



#ifndef UPARSERBEDGRAPH_INCLUDED
#define UPARSERBEDGRAPH_INCLUDED

#include "uParserBase.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>
#include <boost/xpressive/xpressive.hpp>
//#include "uParserFactory.h"
namespace NGS
{

class uParserBedGraph : public uParserBase
{
public :
    uParserBedGraph();
    ~uParserBedGraph();

    virtual void init(const std::string& filename, bool header = true);
    virtual void init(std::istream* stream, bool header = true);
    uToken getNextEntry();
     static uParserBase * Create() { return new uParserBedGraph(); }
private:
    uToken _getTokenFromBedGraphString(const std::string & line);
    bool _parseHeader();
//    static DerivedParserRegister<uParserBedGraph> reg;
    bool m_headerFound=false;
	const std::string s_bedGraphHeader="track type=bedGraph";
    std::vector<std::string> m_tokens;
};

} // End of namespace NGS

#endif // UPARSERBEDGRAPH_INCLUDED
