// ***************************************************************************
// uParser.cpp (c) 2013
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
// **************************************************************************

#include "../../uGeneException.h"
#include "uParser.h"
#include "uParserBase.h"
//#include "uParserBed.h"
#include "uParserFactory.h"
namespace NGS
{

/** \brief Filename default constructor.
 * \param const std::string& filename: the name of the file to parse.
 * \param const std::string& type: the type of file (i.e.: "BED")
 * \param bool header: true if there is a header (value at false by default).
 */
uParser::uParser(const std::string& filename, const std::string & type, bool header)
{

    try
    {
         m_pParserBase=uParserBaseFactory::GetFact()->createInstance(type);
        m_pParserBase->init(filename, header);
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Stream default constructor.
 * \param std::iostream* stream: the stream to parse.
 * \param const std::string& type: the type of stream (i.e.: "BED")
 * \param bool header: true if there is a header (value at false by default).
 */
uParser::uParser(std::istream* stream, const std::string & type, bool header)
{


    m_pParserBase= uParserBaseFactory::GetFact()->createInstance(type);
    m_pParserBase->init(stream, header);

}

/** \brief Filename custom constructor.
 * \param const std::string& filename: the name of the file to parse.
 * \param const std::vector<std::string>& fieldsNames: The name of every column in the file.
 * \param bool header: true if there is a header (value at false by default).
 */
uParser::uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter)
{

    m_pParserBase=uParserBaseFactory::GetFact()->createInstance("CUSTOM");
    m_pParserBase->init(filename, fieldsNames, delimiter);

}

/** \brief Stream custom constructor.
 * \param std::iostream* stream: The stream to parse.
 * \param const std::vector<std::string>& fieldsNames: The name of every column in the file.
 * \param bool header: true if there is a header (value at false by default).
 */
uParser::uParser(std::istream* stream, const std::vector<std::string>& fieldsNames, char delimiter)
{

    //uParserBaseFactory myFact;
     m_pParserBase=uParserBaseFactory::GetFact()->createInstance("CUSTOM");
   // m_pParserBase=myFact.createInstance("CUSTOM");
    m_pParserBase->init(stream, fieldsNames, delimiter);

}

/** \brief Destructor.
 */
uParser::~uParser() {}

/** \brief Check if we are at the end of the file (or of the stream)
 * \return true if we are at the end of the file, otherwise return false.
 */
bool uParser::eof() const
{
    return m_pParserBase->eof();
}

/** \brief Create a token from current point in the file (or the stream).
 * \return a uToken objet containing the infos for the next entry.
 */
uToken uParser::getNextEntry()
{
    return  m_pParserBase->getNextEntry();
}

}
