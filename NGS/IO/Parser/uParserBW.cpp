
// ***************************************************************************
// uParserBW.cpp (c) 2015
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
// **************************************************************************

#include "uParserBW.h"
#include "../../utility/utility.h"
#include <sstream>
//#include "uParserFactory.h"
namespace NGS
{
/** \brief Default empty constructor
 */
uParserBW::uParserBW()
:uParserBase()
{

}
/** \brief Default empty destructor
 */
uParserBW::~uParserBW()
{
}

/** \brief Called when created, loads stream
 *
 * \param filename const std::string& Path to file
 * \return void
 *
 */
void uParserBW::init(const std::string& filename,bool header )
{

    uParserBase::init(filename);
    m_fstream.open(filename);
    m_BBReader.open(filename,m_fstream);
    m_BBIterator = m_BBReader.getBigWigIterator();
}

/** \brief
 *
 * \param stream std::istream*
 * \param header bool leave at false, ignored if set
 * \return void
 *
 */
void uParserBW::init(std::istream* stream, bool header )
{
    throw uParser_exception_base()<<string_error("Calling invalid BWParser version");

   // uParserBase::init(stream, header);
   // _parseHeader();

}



/** \brief Parses next line in Sam file with xpressive regex
 *
 * \return uToken Token created from Sam
 *
 */
uToken uParserBW::getNextEntry()
{
    try
    {
        uToken ourToken;
        ourToken._setParamNoValidate(token_param::CHR,m_BBIterator->getChromosome());
        ourToken._setParamNoValidate(token_param::START_POS,utility::to_string(m_BBIterator->getStartBase()));
        ourToken._setParamNoValidate(token_param::END_POS,utility::to_string(m_BBIterator->getEndBase()));
        ourToken._setParamNoValidate(token_param::SCORE,utility::to_string(m_BBIterator->getWigValue()));

        ++m_BBIterator;

		return ourToken;
		//return uToken(token_infos,false,false);
    }
catch(invalid_uToken_throw& e)
    {
        throw e;
    }
}

bool uParserBW::eof(){

    if ( m_BBIterator.isEnd() ){
         return true;
    }
    return false;
}


//DerivedParserRegister<uParserBW> uParserBW::reg("SAM");
}
