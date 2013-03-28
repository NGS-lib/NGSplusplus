// ***************************************************************************
// uParser.h (c) 2013
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



#ifndef IPARSER_H_INCLUDED
#define IPARSER_H_INCLUDED

#include <map>
#include <string>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "../uToken.h"
#include "../../uGeneException.h"
#include  "uParserBase.h"
#include "../uHeader.h"
namespace NGS {
class uHeader;
class uParser{

public :

    uParser(const std::string& filename, const std::string & type, bool header = false);
    uParser(std::istream* stream, const std::string & type, bool header = false);
    uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    uParser(std::istream* stream, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    ~uParser();

    bool eof() const ;
    uToken getNextEntry();
    std::string getPreviousRaw(){return m_pParserBase->getPreviousRaw();};
    /** \brief Get an unformated version of header (i.e.: a single string containing the whole header)
    */
    std::string getUnformatedHeader() const { return m_pParserBase->getUnformatedHeader(); }
    /** \brief Get a specific data from header.
      */
    std::string getHeaderParam(header_param name) const { return m_pParserBase->getHeaderParam(name); }
    std::vector<std::string> getHeaderParamVector(header_param name) const { return m_pParserBase->getHeaderParamVector(name); }
    /** \brief Check if there is a value associated with a given param.
      * \param header_param& name: name of the param to check.
      */
    bool isHeaderParamSet(const header_param& name) const { return m_pParserBase->isHeaderParamSet(name); }
private:
    std::shared_ptr<uParserBase> m_pParserBase=nullptr;

 //   static uParserBaseFactory myFact;

};
}
#endif
