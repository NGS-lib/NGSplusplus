// ***************************************************************************
// uWriter.h (c) 2013
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



#ifndef UWRITER_H_INCLUDED
#define UWRITER_H_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include <memory>
//#include "../../NGS++.h"
#include "uWriterBase.h"


namespace NGS
{


class uWriter
{
public:
    uWriter(const std::string& filename, const std::string& type);
    uWriter(std::ostream* os, const std::string& type);
    uWriter(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    uWriter(std::ostream* os, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    void writeToken(const uToken& token);
    void printString(const std::string& str);

    void addToHeader(header_param param,std::string value);
    void writeHeader();

private:
    std::shared_ptr<uWriterBase> m_pWriterBase = nullptr;
}; // End of class uWriter

} // End of namespace NGS
#endif // UWRITER_H_INCLUDED
