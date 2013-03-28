// ***************************************************************************
// uParserBase.h (c) 2013
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



#ifndef IUPARSERBASE_H_INCLUDED
#define IUPARSERBASE_H_INCLUDED

#include <map>
#include <string>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "../uToken.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <boost/xpressive/xpressive.hpp>
namespace NGS
{

namespace PDEF{
    const char UCSCCOMMENT='#';
    const std::string UCSCBROWSER="browser";

    inline static bool isUCSCComment(const std::string & pLine){
        if (!pLine.size())
            return false;
        if (pLine.at(0)==(PDEF::UCSCCOMMENT))
            return true;
        else
            return false;
    }
    inline static bool isUCSCBrowser(const std::string & pLine){
     if (pLine.substr(0,7)==PDEF::UCSCBROWSER)
        {
            return true;
        }
    return false;
    }
   inline  static bool isUCSCIgnore(const std::string & pLine)
    {
        if (isUCSCComment(pLine) || isUCSCBrowser(pLine))
            return true;
        else
            return false;
    }
}

class uParserBase
{
    public :
    uParserBase();
    virtual ~uParserBase();

    virtual void init(const std::string& filename, bool header = false);
    virtual void init(std::istream* stream, bool header = false);
    virtual void init(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter = '\t') { };
    virtual void init(std::istream* stream, const std::vector<std::string>& fieldsNames, char delimiter = '\t') { };
    uParserBase& operator=(const uParserBase& copyFrom) = delete;
    uParserBase(const uParserBase&) = delete;
    /** \brief Check if input data is at end of file.
      */
    virtual bool eof();
    virtual uToken getNextEntry()=0;
    std::string getPreviousRaw(){return m_rawString;};
    /** \brief Get an unformated version of header (i.e.: a single string containing the whole header)
    */
    std::string getUnformatedHeader() const { return m_headerData.getUnformatedHeader(); }
    /** \brief Get a specific data from header.
      */
    std::string getHeaderParam(header_param name) const { return m_headerData.getParam(name); } ;
    std::vector<std::string> getHeaderParamVector(header_param name) const{return m_headerData.getParamVector(name); } ;
    /** \brief Check if there is a value associated with a given param.
      * \param header_param& name: name of the param to check.
      */
    bool isHeaderParamSet(const header_param& name) const { return m_headerData.isParamSet(name); }

protected:
    std::stringstream m_hBuffer;
    bool m_dynamicStream = false;
    std::istream* m_pIostream=nullptr;
    uHeader m_headerData;
    std::string m_rawString="";

};

}
#endif // IuParserBase_H_INCLUDED
