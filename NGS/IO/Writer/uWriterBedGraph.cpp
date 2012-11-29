#include "uWriterBedGraph.h"
namespace NGS
{
/** \brief Print the values of a token in BEDGRAPH format in current file
  * \param const uToken& token: the token to print.
  */

void uWriterBedGraph::writeToken(const uToken& token)
{
    try
    {
        /**< Bedgrahp is chr, start, end, score  */
        /**< If no value available, score is 0, thought kind of makes the bedgraph pointless */
        std::string chr="*";
        std::string score="0";
        const std::string TAB = "\t";
        if (token.isParamSet(token_param::CHR))
            chr=token.getParam(token_param::CHR);
        if (token.isParamSet(token_param::SCORE))
            score=token.getParam(token_param::SCORE);
        *m_pOstream << chr<<TAB<< token.getParam(token_param::START_POS)<<TAB << token.getParam(token_param::END_POS)<<TAB <<score <<std::endl;
    }
    catch(param_not_found& e)
    {
        throw e;
    }
}

/** \brief Write basic, theoretically mandatory, bedgraph track type header
 * \return void
 *
 */
void uWriterBedGraph::writeHeader(){

    *m_pOstream<<m_bedGraphHeader<<std::endl;
}


DerivedRegister<uWriterBedGraph> uWriterBedGraph::reg("BEDGRAPH");
} // End of namespace NGS
