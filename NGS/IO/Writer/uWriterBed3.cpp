#include "uWriterBed3.h"

namespace NGS
{

/** \brief Print the values of a token in Bed format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterBed3::writeToken(const uToken& token)
{
    try
    {
        std::string chr="*";
        if (token.isParamSet(token_param::CHR))
        {
             chr = token.getParam(token_param::CHR);
        }
        std::string start_pos = token.getParam(token_param::START_POS);
        std::string end_pos = token.getParam(token_param::END_POS);
        *m_pOstream << chr << '\t' << start_pos << '\t' << end_pos;
    }
    catch(param_not_found& e)
    {
        throw e;
    }
    *m_pOstream << std::endl;
}

//DerivedRegister<uWriterBed3> uWriterBed3::reg("BED3");
} // End of namespace NGS
