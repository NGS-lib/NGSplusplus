#include "uParserBAM.h"
#include "../utility/utility.h"
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
    _parseHeader();

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
        ourToken._setParamNoValidate(token_param::FLAGS, std::to_string(m_BufferAlignement.AlignmentFlag) );
        ourToken._setParamNoValidate(token_param::CHR,  m_BamReader.GetReferenceData().at(m_BufferAlignement.RefID).RefName);
        ourToken._setParamNoValidate(token_param::START_POS, std::to_string(m_BufferAlignement.Position));
        ourToken._setParamNoValidate(token_param::MAP_SCORE, std::to_string(m_BufferAlignement.MapQuality));
        ourToken._setParamNoValidate(token_param::TEMPLATE_LENGHT, std::to_string(m_BufferAlignement.InsertSize));
        ourToken._setParamNoValidate(token_param::SEQUENCE,   m_BufferAlignement.QueryBases );
        ourToken._setParamNoValidate(token_param::PHRED_SCORE, m_BufferAlignement.Qualities);
        std::string cigar;
        for(CigarOp & cigarItem:  m_BufferAlignement.CigarData)
        {
           cigar+= ( cigarItem.Type+std::to_string(cigarItem.Length));
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
        return true;

    if ( m_BamReader.GetNextAlignment(m_BufferAlignement) )
        m_IsBuffer=true;
        return true;
    return false;
}


//DerivedParserRegister<uParserBAM> uParserBAM::reg("SAM");
}
