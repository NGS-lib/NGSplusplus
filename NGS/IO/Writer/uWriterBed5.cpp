#include "uWriterBed5.h"

namespace NGS
{

/** \brief Print the values of a token in Bed format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterBed5::writeToken(const uToken& token)
{
    try
    {
        uWriterBed::writeToken(token);
        /**< If score is not set, we replace it with the value "." */
        std::string score=".";
        if (token.isParamSet(token_param::SCORE))
        {
            score = token.getParam(token_param::SCORE);
        }
        *m_pOstream << '\t' << score;
    }
    catch(param_not_found& e)
    {
        throw e;
    }
    *m_pOstream << std::endl;
}
} // End of namespace NGS
