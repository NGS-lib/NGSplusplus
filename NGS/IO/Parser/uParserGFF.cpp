#include "uParserGFF.h"

namespace NGS {

/* Note, there are ambiguities in the GFF format. The spec at http://www.sanger.ac.uk/resources/software/gff/spec.html
describes the first column as

"The name of the sequence. Having an explicit sequence name allows a feature file to be prepared for a data set of multiple sequences.
Normally the seqname will be the identifier of the sequence in an accompanying fasta format file. An alternative is that <seqname>
is the identifier for a sequence in a public database, such as an EMBL/Genbank/DDBJ accession number. Which is the case, and which
file or database to use, should be explained in accompanying information. "


However, practically, use tend to set column one as scaffp;d ID ( ex: chromosome ). This same langage is used for GFF3

As such, we parser SEQ_NAME from first column and do not set CHR

Note that the GFF3 parser uses the first column as CHR
 */



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
    /**< GFF regex */
     GFFRegex = sregex::compile(GFFregString) ;
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

uToken uParserGFF::_getTokenFromGFFString(const std::string & line)
{
    /**< This would be more efficient with a static regex, but preserving Perl syntax makes it "easier" to read */
    smatch what;
    if( regex_match( line, what, GFFRegex ) )
    {
        /**< Preset according to GFF2 format */
        uToken ourToken;
        /**< According to GFF specification, the first value is "SEQNAME". Howeverm, practical use for most people is using it as chrom.  */
        /**< As such, the assignation of what[1] is subject to change */
        if ( what[1]!=".")
        ourToken._setParamNoValidate(token_param::SEQ_NAME, what[1]);
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
        if ( what[8]!=".")
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


//DerivedParserRegister<uParserGFF> uParserGFF::reg("GFF");
} // End of namespace NGS
