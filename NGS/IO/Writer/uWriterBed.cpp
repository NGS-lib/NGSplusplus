#include "uWriterBed.h"

namespace NGS
{

/** \brief Print the values of a token in Bed format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterBed::writeToken(const uToken& token)
{
    try
    {
        std::string chr="*";
        if (token.isParamSet(token_param::CHR))
             chr = token.getParam(token_param::CHR);
        std::string start_pos = token.getParam(token_param::START_POS);
        std::string end_pos = token.getParam(token_param::END_POS);
        /**< If seq_name is not set, we replace it with the value "." */
        std::string seq_name=".";
        if (token.isParamSet(token_param::SEQ_NAME))
        {
            seq_name = token.getParam(token_param::SEQ_NAME);
        }
        *m_pOstream << chr << '\t' << start_pos << '\t' << end_pos << '\t' << seq_name;
    }
    catch(param_not_found& e)
    {
        throw e;
    }
}

} // End of namespace NGS
