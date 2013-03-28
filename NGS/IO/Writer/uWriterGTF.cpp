// ***************************************************************************
// uWriterGTF.cpp (c) 2013
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


#include "uWriterGTF.h"
namespace NGS
{
/** \brief Print the values of a token in GFF format in current file
  * \param const uToken& token: the token to print.
  */


/*
Practically for now, the onlye difference between this and GFF is that we assign column 1 to CHR rather then SEQ_NAME
 */
void uWriterGTF::writeToken(const uToken& token)
{
    try
    {
        /**< If score or strand are not set, we replace them with the value "." */
        const std::string TAB="\t";
        std::string seqname=".";
        std::string feature=".";
        std::string score=".";
        std::string strand=".";
        std::string start;
		std::string end;
		std::string source=".";
		std::string phase=".";
		std::string extra, group_ID="", group_transcript="";

		/**< No default values */
		start=token.getParam(token_param::START_POS);
		end=token.getParam(token_param::END_POS);

        if (token.isParamSet(token_param::SCORE))
       	 { score = token.getParam(token_param::SCORE); }

		if(token.isParamSet(token_param::CHR))
       		 { seqname = token.getParam(token_param::CHR); }
        else
		if (token.isParamSet(token_param::SEQ_NAME))
            { seqname = token.getParam(token_param::SEQ_NAME); }

		if(token.isParamSet(token_param::FEATURE_TYPE))
			{ feature = token.getParam(token_param::FEATURE_TYPE); }
        /**< If score is not set, but strand is, we replace score by "." */
        if (token.isParamSet(token_param::STRAND))
        {  strand = token.getParam(token_param::STRAND);   }

		if (token.isParamSet(token_param::SOURCE))
        {  source = token.getParam(token_param::SOURCE);   }

		if(token.isParamSet(token_param::PHASE))
			{ phase = token.getParam(token_param::PHASE); }

        if(token.isParamSet(token_param::PHASE))
			{ phase = token.getParam(token_param::PHASE); }

        if(token.isParamSet(token_param::GROUP_ID))
			{ group_ID = token.getParam(token_param::GROUP_ID); }

        if(token.isParamSet(token_param::GROUP_TRANSCRIPT))
			{ group_transcript = token.getParam(token_param::GROUP_TRANSCRIPT); }

		if(token.isParamSet(token_param::EXTRA))
       		 { extra = token.getParam(token_param::EXTRA); }

        *m_pOstream << seqname<<TAB<< source<<TAB << feature<<TAB << start<<TAB<<end <<TAB << score<<TAB<<strand <<TAB<< phase<<TAB<<"gene_id \""<<group_ID<<"\"; "<<"transcript_id \""<<group_transcript<<"\";\n";
    }
    catch(param_not_found& e)
    {
        throw e;
    }
}
//DerivedRegister<uWriterGTF> uWriterGTF::reg("GTF");
} // End of namespace NGS
