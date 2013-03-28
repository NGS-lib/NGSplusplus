// ***************************************************************************
// uParserBase.cpp (c) 2013
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


#include "uParserBase.h"
#include "uParserFactory.h"
namespace NGS
{

/** \brief Destructor.
 */
uParserBase::~uParserBase()
{
    if (m_dynamicStream == true)
    {
        delete m_pIostream;
    }
    m_pIostream = NULL;
}

/** \brief Initialize the uParserBase object (open file and parse header).
 * \param const std::string& filename: name of the bed file to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserBase::init(const std::string& filename, bool header)
{
    std::ifstream* ifs = new std::ifstream(filename.c_str(), std::ifstream::in);
    if (!ifs->is_open())
    {
        std::string error = "Error opening file: " + filename;
        throw  std::runtime_error(error.c_str());
    }
    else
    {
        m_pIostream = ifs;
        m_dynamicStream = true;
    }
}

/** \brief Initialize the uParserBase object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserBase::init(std::istream* stream, bool header)
{
    m_pIostream = stream;
    m_dynamicStream = false;
}


/** \brief Default constructor.
 */
uParserBase::uParserBase(){}

/** \brief Check if we are at the end of the file/stream including if the last line concludes with a NL
 * \return true is we are at the end, otherwise false.
 */
bool uParserBase::eof()
{
    if (m_pIostream->peek() == EOF)
        return true;
    if (m_pIostream->peek()=='\n')
    {
       m_pIostream->get();
       if (m_pIostream->peek() ==EOF)
           return true;
       else
       {
        m_pIostream->unget();
        return false;
       }
    }
    return m_pIostream->peek() == EOF;

}

std::map<std::string, std::function<uParserBase*()> > *uParserBaseFactory::mapItem;
}
