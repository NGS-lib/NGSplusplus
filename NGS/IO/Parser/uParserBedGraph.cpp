#include "../../utility/utility.h"
#include "uParserBedGraph.h"
namespace NGS
{

/** \brief Default constructor (not used directly).
 */
uParserBedGraph::uParserBedGraph(): uParserBase()
{
}
/** \brief Destructor.
 */
uParserBedGraph::~uParserBedGraph()
{
}
/** \brief Initialize the uParserBedGraph object (open file and parse header).
 * \param const std::string& filename: name of the bed file to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserBedGraph::init(const std::string& filename, bool header)
{
    uParserBase::init(filename, header);
    _parseHeader();
}

/** \brief Initialize the uParserBedGraph object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserBedGraph::init(std::istream* stream, bool header)
{
    uParserBase::init(stream, header);
    _parseHeader();

}

/** \brief Produce a token with next entry in the file/stream.
 * \return uToken containing the infos of the next entry.
 */
uToken uParserBedGraph::getNextEntry()
{
    std::string strLine;
    /**< Do we need to parse our buffer */
    if (m_hBuffer.eof()==false)
        std::getline(m_hBuffer, strLine);
        m_rawString=strLine;
        /**< If we stored something, otherwise, getline, otherwise fail */
        if ( (strLine.size()) || ( std::getline(*m_pIostream, strLine)) )
        {
            /**< We start by fetching the infos from the line */
            return _getTokenFromBedGraphString(strLine);
        }
        else
        {
#ifdef DEBUG
            std::cerr << "Reached end of file." << std::endl;
#endif
            end_of_file_throw e;
            e << string_error("Reached end of file. BedGraph header object");
            throw e;
        }
    std::cerr <<"Fatal error in getNextEntry() from BedGraph parser, should not reach here." <<std::endl;
    abort();
}

uToken uParserBedGraph::_getTokenFromBedGraphString(const std::string & line)
{
    ///**< A Bedgraph format is <string> <int> <int> <float>, simple enough that we skip regex */
    try
    {
        utility::GetTokens(m_tokens,line);
        /**< Validate */
        if (!(utility::is_posnumber(m_tokens.at(1) ))&&(utility::is_posnumber(m_tokens.at(2)))&&(m_tokens.size()==4))
            throw;
        m_tokens.at(3).c_str();;
        uToken ourToken;
        ourToken._setParamNoValidate(token_param::CHR,m_tokens.at(0));
        ourToken._setParamNoValidate(token_param::START_POS,m_tokens.at(1));
        ourToken._setParamNoValidate(token_param::END_POS,m_tokens.at(2));
        ourToken._setParamNoValidate(token_param::SCORE,m_tokens.at(3));
        return ourToken;
      //  token_infos << "CHR\t" << ss << "\n";
      //  token_infos << "START_POS\t" << ss << "\n";
      //  token_infos << "END_POS\t" << ss << "\n";
     //   token_infos << "SCORE\t" << ss << "\n";
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Failed while parsing bedgraph line" << std::endl;
#endif
        uParser_invalid_BedGraph_line e;
        e << string_error("Failed while parsing the following bedgraph line:\n"+line);
        throw e;
    }

}

/** \brief Parse the (theoretically) mandatory bedgraph header, report if you found it and skip line
 * \return void
 */
bool uParserBedGraph::_parseHeader()
{
    /**< Note, we cannot validate there is a header simply by peaking.  */
    /**< So we load the line and if it is not a header, store it in our buffer object for parsing */
    /**< BedGraph header information is very specified to UCSC, so we do not store it. */
    std::string strLine;
    if (std::getline(*m_pIostream, strLine))
    {
        /**< If we could not find track type=bedGraph, store, otherwise discard */
        if (strLine.find(s_bedGraphHeader)==std::string::npos){
            m_hBuffer << strLine;
            return false;
        }
        return true;
    }
    else
    {
#ifdef DEBUG
        std::cerr << "Reached end of file." << std::endl;
#endif
        end_of_file_throw e;
        e << string_error("Reached end of file when parsing bedGraph header.");
        throw e;
    }
}

//DerivedParserRegister<uParserBedGraph> uParserBedGraph::reg("BEDGRAPH");
} // End of namespace NGS
