#include "uParserGTF.h"

namespace NGS {


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

    /**< GTF regex */
     sregex GTFRegex = sregex::compile(GTFregString) ;
}

/** \brief Initialize the uParserGTF object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserGTF::init(std::istream* stream, bool header)
{
    uParserBase::init(stream, header);
    /**< GTF regex */
     sregex GTFRegex = sregex::compile(GTFregString) ;
}

/** \brief Produce a token with next entry in the file/stream.
 * \return uToken containing the infos of the next entry.
 */
uToken uParserGTF::getNextEntry()
{
  std::string strLine;
   // char line[4096];
    //if (m_pIostream->getline(strLine))
    if (std::getline(*m_pIostream, strLine))
    {
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
       // token_infos << "SEQ_NAME\t" << what[1] << "\n";
        if ( what[2]!=".")
        ourToken._setParamNoValidate(token_param::SOURCE, what[2]);
      //  token_infos << "SOURCE\t" << what[2] << "\n";
       if ( what[3]!=".")
        ourToken._setParamNoValidate(token_param::FEATURE_NAME, what[3]);
       // token_infos << "FEATURE_NAME\t" << what[3] << "\n";
        ourToken._setParamNoValidate(token_param::START_POS, what[4]);
       // token_infos << "START_POS\t" << what[4] << "\n";
       ourToken._setParamNoValidate(token_param::END_POS, what[5]);
       // token_infos << "END_POS\t" <<  what[5] << "\n";
        if ( what[6]!=".")
             ourToken._setParamNoValidate(token_param::SCORE, what[6]);
            //token_infos << "SCORE\t" << what[6] << "\n";
        /**< GFF considered a '.' to mean no info or not relevant. We simply do not stock it */
        if ( what[7]!=".")
            ourToken._setParamNoValidate(token_param::STRAND, what[7]);
           // token_infos << "STRAND\t" << what[7] << "\n";
             ourToken._setParamNoValidate(token_param::PHASE, what[8]);
        //token_infos << "PHASE\t" << what[8] << "\n";
        if (what[9].matched)
             ourToken._setParamNoValidate(token_param::EXTRA, what[9]);
            //token_infos << "EXTRA\t" << what[9] << "\n";
        return ourToken;
    }
    else{
        throw uParser_invalid_GFF_line()<<string_error("GFF line, failling validation. Line is:\n"+line);
    }

}

//DerivedParserRegister<uParserGTF> uParserGTF::reg("GTF");
} // End of namespace NGS
