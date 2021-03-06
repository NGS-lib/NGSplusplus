// ***************************************************************************
// uParserGenePred.cpp (c) 2013
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


#include "uParserGenePred.h"
namespace NGS
{

using namespace boost::xpressive;
/** \brief Default constructor (not used directly).
 */

uParserGenePred::uParserGenePred(): uParserBase()
{
}

/** \brief Destructor.
 */
uParserGenePred::~uParserGenePred()
{
}

/** \brief Initialize the uParserGenePred object (open file and parse header).
 * \param const std::string& filename: name of the bed file to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserGenePred::init(const std::string& filename, bool header)
{
    uParserBase::init(filename, header);

    /**< GFF regex */
    UCSCGenePredRegexInit = sregex::compile(GenePredregStringPrior) ;
    _parseHeader();
}

/** \brief Initialize the uParserGenePred object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserGenePred::init(std::istream* stream, bool header)
{
    uParserBase::init(stream, header);
    /**< GFF regex */
    UCSCGenePredRegexInit = sregex::compile(GenePredregStringPrior) ;
    _parseHeader();
}

 void uParserGenePred::_parseHeader()
 {
    /**< skip header lines */
    std::string strLine;
    while (!m_pIostream->eof() &&  m_pIostream->peek()=='#')
        {std::getline(*m_pIostream, strLine); }
 }


/** \brief Produce a token with next entry in the file/stream.
 * \return uToken containing the infos of the next entry.
 */
uToken uParserGenePred::getNextEntry()
{
    std::string strLine;

    if  ( std::getline(*m_pIostream, strLine) && ( strLine.size() || m_pIostream->eof()==false) )
    {
        m_rawString=strLine;
        return _getTokenFromGenePredString(strLine);
    }
    else
    {
#ifdef DEBUG
        std::cerr << "Reached end of file." << std::endl;
#endif
        end_of_file_throw e;
        e << string_error("Reached end of file.");
        throw e;
    }
    std::cerr <<"Fatal error in getNextEntry() from GFF parser, should not reach here." <<std::endl;
    abort();
}

uToken uParserGenePred::_getTokenFromGenePredString(const std::string & line)
{
    /**< This would be more efficient with a static regex, but preserving Perl syntax makes it "easier" to read */
    smatch what;
    if( regex_match( line, what, UCSCGenePredRegexInit ) )
    {
        /**< GenePred normal Format*/

        /**< Seize number of exons */
        int exonCount=utility::stoi(what[8]);
        std::string completeRegex=GenePredregStringPart1;
        /**< Make complete regex */
        if (exonCount){
            completeRegex+="\t";
            for (int i=0; i<exonCount; i++)
            {
                completeRegex+=GenePredRegExon;
            }
            completeRegex+="\t";
            for (int i=0; i<exonCount; i++)
            {
                completeRegex+=GenePredRegExon;
            }

        }
        /**< To make compatible with any number of genePred extensions */
        completeRegex+="(.*)";
        UCSCGenePredRegex = sregex::compile(completeRegex) ;
        if( regex_match( line, what, UCSCGenePredRegex ) )
        {
            uToken ourToken;

            ourToken._setParam(token_param::GROUP_ID, what[1]);
            ourToken._setParam(token_param::CHR, what[2]);
            ourToken._setParam(token_param::STRAND, what[3]);
            ourToken._setParam(token_param::START_POS, what[4]);
            ourToken._setParam(token_param::END_POS, what[5]);
            ourToken._setParam(token_param::FEATURE_TYPE, "GENE");
            ourToken._setParam(token_param::START_POS, what[6]);
            ourToken._setParam(token_param::END_POS, what[7]);
            ourToken._setParam(token_param::FEATURE_TYPE, "CODING");
            for (int i=0; i<exonCount; i++)
            {
                ourToken._setParam(token_param::START_POS, what[i+9]);
                ourToken._setParam(token_param::END_POS, what[i+9+exonCount]);
                ourToken._setParam(token_param::FEATURE_TYPE, "EXON");
            }
            return ourToken;
        }
        else
        {
            throw uParser_invalid_GenePred_line()<<string_error("GFF line, failling validation. Line is:\n"+line);

        }
    }
    else
    {
        throw uParser_invalid_GFF_line()<<string_error("GFF line, failling validation. Line is:\n"+line);

    }
}

} // End of namespace NGS
