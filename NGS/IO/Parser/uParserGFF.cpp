#include "uParserGFF.h"

namespace NGS {


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

    m_delimiter = '\t';

    /**< GFF regex */
     sregex GFFRegex = sregex::compile(GFFregString) ;
}

/** \brief Initialize the uParserGFF object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserGFF::init(std::istream* stream, bool header)
{
    uParserBase::init(stream, header);
    /**< GFF regex */
     sregex GFFRegex = sregex::compile(GFFregString) ;
}

/** \brief Produce a token with next entry in the file/stream.
 * \return uToken containing the infos of the next entry.
 */
uToken uParserGFF::getNextEntry()
{
    std::string strLine;
   // char line[4096];
    //if (m_pIostream->getline(strLine))
    if (std::getline(*m_pIostream, strLine))
    {
        /**< We start by fetching the infos from the line */
        std::stringstream token_infos;
        _getTokenInfoFromGFFString(strLine, token_infos);
        /**< We try to create a token with the infos that were fetched from the line */
        try
        {
            return uToken(token_infos);
        }
        catch(invalid_uToken_throw& e)
        {
            throw e;
        }
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

void uParserGFF::_getTokenInfoFromGFFString(std::string line, std::stringstream& token_infos)
{
    /**< Matchin GFF 2 format, to following */
    /**< string string string int int float +/-/. int  */
    /**< So minor post-regex validation needed */
    /**< This would be more efficient with a static regex, but preserving Perl syntax makes it "easier" to read */
    smatch what;
     boost::xpressive::sregex GFFRegex2=sregex::compile(GFFregString) ;


	//std::string test("SEQ1\tEMBL\tsplice5\t172\t173\t.\t+\t.");

     sregex GFFRegex3 = sregex::compile("^([\\w_-]+)\t([\\w_-]+)\t([\\w_-]+)\t(\\d+)\t(\\d+)\t([-+]?[0-9]*\\.?[0-9]+|.)\t(\\+|\\-|\\.)\t([012\\.])(?:\t(.+))?") ;
    if( regex_match( line, what, GFFRegex3 ) )
    {
        /**< Preset according to GFF2 format */

        /**< According to GFF specification, the first value is "SEQNAME". Howeverm, practical use for most people is using it as chrom.
        /**< As such, the assignation of what[1] is subject to change */

        token_infos << "SEQ_NAME\t" << what[1] << "\n";
        token_infos << "SOURCE\t" << what[2] << "\n";
        token_infos << "FEATURE_NAME\t" << what[3] << "\n";
        token_infos << "START_POS\t" << what[4] << "\n";
        token_infos << "END_POS\t" <<  what[5] << "\n";
        if ( what[6]!=".")
            token_infos << "SCORE\t" << what[6] << "\n";
        /**< GFF considered a '.' to mean no info or not relevant. We simply do not stock it */
        if ( what[7]!=".")
            token_infos << "STRAND\t" << what[7] << "\n";
        token_infos << "PHASE\t" << what[8] << "\n";
        if (what[9].matched)
            token_infos << "EXTRA\t" << what[9] << "\n";
    }
    else{
        throw uParser_invalid_GFF_line()<<string_error("GFF line, failling validation. Line is:\n"+line);
    }

}


DerivedParserRegister<uParserGFF> uParserGFF::reg("GFF");
} // End of namespace NGS
