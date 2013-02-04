#include "uWriterBed4.h"

namespace NGS
{

/** \brief Print the values of a token in Bed format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterBed4::writeToken(const uToken& token)
{
    try
    {
        uWriterBed::writeToken(token);
    }
    catch(param_not_found& e)
    {
        throw e;
    }
    *m_pOstream << std::endl;
}

//DerivedRegister<uWriterBed4> uWriterBed4::reg("BED4");
} // End of namespace NGS
