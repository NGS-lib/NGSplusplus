// ***************************************************************************
// UWriterBedGraph.cpp (c) 2013
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


//DerivedRegister<uWriterBedGraph> uWriterBedGraph::reg("BEDGRAPH");
} // End of namespace NGS
