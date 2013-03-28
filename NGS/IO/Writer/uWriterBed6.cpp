// ***************************************************************************
// uWriterBed6.cpp (c) 2013
// Alexei Nordell-Markovits : Sherbrooke University
// Charles Joly Beauparlant : Laval University
//
//       This file is part of the NGS++ library.
//
//    The NGS++ library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program (lgpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// **************************************************************************


#include "uWriterBed6.h"

namespace NGS
{

/** \brief Print the values of a token in Bed format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterBed6::writeToken(const uToken& token)
{
    try
    {
        uWriterBed::writeToken(token);
        /**< If score or strand are not set, we replace them with the value "." */
        std::string score=".";
        std::string strand=".";
        if (token.isParamSet(token_param::SCORE))
        {
            score = token.getParam(token_param::SCORE);
        }

        /**< If score is not set, but strand is, we replace score by "." */
        if (token.isParamSet(token_param::STRAND))
        {
            strand = token.getParam(token_param::STRAND);
        }

        *m_pOstream << '\t' << score << '\t' << strand;
    }
    catch(param_not_found& e)
    {
        throw e;
    }
    *m_pOstream << std::endl;
}

//DerivedRegister<uWriterBed6> uWriterBed6::reg("BED6");
} // End of namespace NGS
