// ***************************************************************************
// uParserBAM.cpp (c) 2013
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

#include "uParserBAM.h"
#include "../../utility/utility.h"
#include <sstream>
//#include "uParserFactory.h"
namespace NGS
{
using namespace BamTools;
using namespace boost::xpressive;
/** \brief Default empty constructor
 */
uParserBAM::uParserBAM():uParserBase()
{

}
/** \brief Default empty destructor
 */
uParserBAM::~uParserBAM()
{
}

/** \brief Called when created, loads stream and parses BAM file header, loading mandatory info
 *
 * \param filename const std::string& Path to file
 * \param header bool leave at false, ignored if set
 * \return void
 *
 */
void uParserBAM::init(const std::string& filename, bool header )
{
   // uParserBase::init(filename, header);
    m_BamReader.Open(filename);
      //  throw parser_exception
  //  _parseHeader();

//            m_headerData._addToParam(header_param::CHR,chrom);
 //           m_headerData._addToParam(header_param::CHR_SIZE,utility::to_string(refSeqlenght));
}
/** \brief Called when created, reads from stream and parses BAM file header, loading mandatory info
 *
 * \param stream std::istream*
 * \param header bool leave at false, ignored if set
 * \return void
 *
 */
void uParserBAM::init(std::istream* stream, bool header )
{
    throw uParser_exception_base()<<string_error("Calling invalid bamParser version");

   // uParserBase::init(stream, header);
   // _parseHeader();

}


/** \brief Parses next line in Sam file with xpressive regex
 *
 * \return uToken Token created from Sam
 *
 */
uToken uParserBAM::getNextEntry()
{
    try
    {
        /**< if no buffer, load data */
        if (!m_IsBuffer){
            if (!(m_BamReader.GetNextAlignment(m_BufferAlignement))){
                throw end_of_file_throw()<<string_error("Reached eof when reading BAM alignement in uParserBAM::getNextEntry()");
            }
        }
        uToken ourToken;
        /**< Make token from the Bam entry */
        ourToken._setParamNoValidate(token_param::SEQ_NAME, m_BufferAlignement.Name);
        ourToken._setParamNoValidate(token_param::START_POS, utility::to_string(m_BufferAlignement.Position+1));
        ourToken._setParamNoValidate(token_param::FLAGS, utility::to_string(m_BufferAlignement.AlignmentFlag) );
        ourToken._setParamNoValidate(token_param::CHR,  m_BamReader.GetReferenceData().at(m_BufferAlignement.RefID).RefName);
        ourToken._setParamNoValidate(token_param::MAP_SCORE, utility::to_string(m_BufferAlignement.MapQuality));
        ourToken._setParamNoValidate(token_param::TEMPLATE_LENGTH, utility::to_string(m_BufferAlignement.InsertSize));
        ourToken._setParamNoValidate(token_param::SEQUENCE,   m_BufferAlignement.QueryBases );
        ourToken._setParamNoValidate(token_param::PHRED_SCORE, m_BufferAlignement.Qualities);
        std::string cigar;
        for(CigarOp & cigarItem:  m_BufferAlignement.CigarData)
        {
           cigar+= ( cigarItem.Type+utility::to_string(cigarItem.Length));
        }
        ourToken._setParamNoValidate(token_param::CIGAR, cigar);

        m_IsBuffer=false;
        return ourToken;
    }
catch(invalid_uToken_throw& e)
{
    throw e;
}
}

/** \brief Parse SAM header, if any. Store in member variables
 *
 * \return void
 *
 */
void uParserBAM::_parseHeader()
{


}

bool uParserBAM::eof(){

    if (m_IsBuffer)
        return false;

    if ( m_BamReader.GetNextAlignment(m_BufferAlignement) ){
        m_IsBuffer=true;
        return false;
    }
    return true;
}


//DerivedParserRegister<uParserBAM> uParserBAM::reg("SAM");
}
