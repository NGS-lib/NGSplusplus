// ***************************************************************************
// uParserGFF.cpp (c) 2013
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


#include "uParserGFF.h"

namespace NGS
{

using namespace boost::xpressive;
/** \brief Default constructor (not used directly).
 */

uParserGFF::uParserGFF(): uParserBase()
{
}

/** \brief Destructor.
 */
uParserGFF::~uParserGFF()
{
}

/** \brief Initialize the uParserGFF object (open file and parse header).
 * \param const std::string& filename: name of the bed file to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserGFF::init(const std::string& filename, bool header)
{
    uParserBase::init(filename, header);
    _parseHeader();
    /**< GFF regex */
    GFFRegex = sregex::compile(GFFregString) ;
}

/** \brief Initialize the uParserGFF object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserGFF::init(std::istream* stream, bool header)
{
    uParserBase::init(stream, header);
    _parseHeader();
    /**< GFF regex */
    GFFRegex = sregex::compile(GFFregString) ;
}

/** \brief Produce a token with next entry in the file/stream.
 * \return uToken containing the infos of the next entry.
 */
uToken uParserGFF::getNextEntry()
{
    std::string strLine;

    if (m_hBuffer.eof()==false)
            std::getline(m_hBuffer, strLine);
//&& ( strLine.size() || m_pIostream->eof()==false)
    if (strLine.size() || std::getline(*m_pIostream, strLine)  )
    {
        m_rawString=strLine;
        return _getTokenFromGFFString(strLine);
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


void uParserGFF::_parseHeader()
{
   /**< While not specific to the GTF format, we will skip header lines that are UCSC browser tags */
    std::string strLine;
    while(std::getline(*m_pIostream, strLine))
    {
        /**< Skip comment and browser lines, if not check if track line */
        if (PDEF::isUCSCIgnore(strLine)==false)
        {
            /**< No longer a commentary, store in buffer */
          //  std::cout <<"Buffer is"<<m_hBuffer<<'\n';
            m_hBuffer << strLine;
            break;
        }
    }
}



uToken uParserGFF::_getTokenFromGFFString(const std::string & line)
{
    /**< This would be more efficient with a static regex, but preserving Perl syntax makes it "easier" to read */
    smatch what;
    if( regex_match( line, what, GFFRegex ) )
    {
        /**< Preset according to GFF version 2  format as defined here */
        /**<  http://www.sanger.ac.uk/resources/software/gff/spec.html */

        uToken ourToken;
        /**< As such, the assignation of what[1] is subject to change */

        /**< GFF considered a '.' to mean no info or not relevant. We simply do not stock it */
        if ( what[1]!=".")
            ourToken._setParamNoValidate(token_param::CHR, what[1]);

        if ( what[2]!=".")
            ourToken._setParamNoValidate(token_param::SOURCE, what[2]);

        if ( what[3]!=".")
            ourToken._setParamNoValidate(token_param::FEATURE_TYPE, what[3]);

        ourToken._setParamNoValidate(token_param::START_POS, what[4]);
        ourToken._setParamNoValidate(token_param::END_POS, what[5]);

        if ( what[6]!=".")
            ourToken._setParamNoValidate(token_param::SCORE, what[6]);

        if ( what[7]!=".")
            ourToken._setParamNoValidate(token_param::STRAND, what[7]);

        if ( what[8]!=".")
            ourToken._setParamNoValidate(token_param::PHASE, what[8]);

        if (what[9].matched)
            ourToken._setParamNoValidate(token_param::GROUP_ID, what[9]);

        return ourToken;
    }
    else
    {
        throw uParser_invalid_GFF_line()<<string_error("GFF line, failling validation. Line is:\n"+line);
    }
}


//DerivedParserRegister<uParserGFF> uParserGFF::reg("GFF");
} // End of namespace NGS
