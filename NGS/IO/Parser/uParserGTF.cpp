#include "uParserGTF.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
namespace NGS
{


using namespace boost::xpressive;
/** \brief Default constructor (not used directly).
 */

uParserGTF::uParserGTF(): uParserBase()
{
}

/** \brief Destructor.
 */
uParserGTF::~uParserGTF()
{
}

/** \brief Initialize the uParserGTF object (open file and parse header).
 * \param const std::string& filename: name of the bed file to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserGTF::init(const std::string& filename, bool header)
{
    uParserBase::init(filename, header);
    _parseHeader();
    /**< GTF regex */
    GTFRegex = sregex::compile(GTFregString) ;
}

/** \brief Initialize the uParserGTF object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserGTF::init(std::istream* stream, bool header)
{
    uParserBase::init(stream, header);
    _parseHeader();
    /**< GTF regex */
    GTFRegex = sregex::compile(GTFregString) ;
}

/** \brief Produce a token with next entry in the file/stream.
 * \return uToken containing the infos of the next entry.
 */
uToken uParserGTF::getNextEntry()
{
    std::string strLine;
    if (m_hBuffer.eof()==false)
        std::getline(m_hBuffer, strLine);

    if  (strLine.size() || std::getline(*m_pIostream, strLine) )
    {
        m_rawString=strLine;
        return _getTokenInfoFromGTFString(strLine);
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

uToken uParserGTF::_getTokenInfoFromGTFString(const std::string& line)
{
    smatch what;
    if( regex_match( line, what, GTFRegex ) )
    {
        /**< Preset according to GTF format */
        uToken ourToken;
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

        /**< GFF considered a '.' to mean no info or not relevant. We simply do not stock it */
        if ( what[7]!=".")
            ourToken._setParamNoValidate(token_param::STRAND, what[7]);

        ourToken._setParamNoValidate(token_param::PHASE, what[8]);

        if (what[9].matched)
            ourToken._setParamNoValidate(token_param::EXTRA, what[9]);

        return ourToken;
    }
    else
    {
        throw uParser_invalid_GFF_line()<<string_error("GFF line, failling validation. Line is:\n"+line);
    }

}

void uParserGTF::_parseHeader()
{
   /**< While not specific to the GTF format, we will skip header lines that are UCSC browser tags */
    std::string strLine;
    while(std::getline(*m_pIostream, strLine))
    {
        /**< Skip comment and browser lines, if not check if track line */
        if (PDEF::isUCSCIgnore(strLine)==false)
        {
            /**< No longer a commentary, store in buffer */
            m_hBuffer << strLine;
            break;
        }
    }
}

} // End of namespace NGS
