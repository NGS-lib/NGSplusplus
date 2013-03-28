// ***************************************************************************
// uParserCustom.h (c) 2013
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



#ifndef UPARSERCUSTOM_H_INCLUDED
#define UPARSERCUSTOM_H_INCLUDED

#include "uParserBed.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>
//#include "uParserFactory.h"
namespace NGS
{

class uParserCustom : public uParserBed
{
public :
    uParserCustom();
    ~uParserCustom();

    virtual void init(const std::string& filename, bool header = false);
    virtual void init(std::istream* stream, bool header = false);
    void init(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    void init(std::istream* stream, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    static uParserBase * Create() { return new uParserCustom(); }
    virtual uToken getNextEntry();

private:
    std::vector<std::string> m_customFieldNames{};
    void _customParserValidateFields(const std::vector<std::string>& fieldsNames) const;
    void _customParserCopyFields(const std::vector<std::string>& fieldsNames);
    void _convertLineToTokenInfosCustom(char* line, std::stringstream& token_infos);
    bool _paramExists(const std::string& name, const std::vector<std::string>& list) const;

};

} // End of namespace NGS
#endif // UPARSERCUSTOM_H_INCLUDED
