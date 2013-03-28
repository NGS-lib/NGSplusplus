// ***************************************************************************
// uWriterBed3.cpp (c) 2013
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
