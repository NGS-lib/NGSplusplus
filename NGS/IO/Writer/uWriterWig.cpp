// ***************************************************************************
// uWriterWig.cpp (c) 2013
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


#include "uWriterWig.h"

namespace NGS {
using namespace std;


/** \brief The write will fetch the type specified when setting the "score" of the wig file.
 *
 * \param type const token_param Value type to write
 * \return void
 *
 */
void uWriterWig::setValueType(const token_param type){

    m_scoreType=type;
}

/** \brief Print the values of a token in SAM format in current file
  * \param const uToken& token: the token to print.
   */
  void uWriterWig::writeToken(const uToken& token) {
 	try {
       if  ((token.isParamSet(token_param::CHR))&&(token.isParamSet(token_param::START_POS))&&(token.isParamSet(token_param::END_POS)) )
           {
                long int ourstart= utility::stoll(token.getParam(token_param::START_POS));
                long int ourend= utility::stoll(token.getParam(token_param::END_POS));
                int ourspan = ( ourend-ourstart+1);
                float density =1;
                if (token.isParamSet(m_scoreType))
                {
                 /**< Should not be using catch as a control flow method, but will to for now */
                 /**< Establish standard for bed conversion */
                 try {
                    density=utility::stof(token.getParam(m_scoreType));
                 }
                    catch(...){
                    density =1;
                    }

                }
                const string step="variableStep chrom=";
                const string span="span=";
                bool writeHeader=false;
                if (m_chr!=token.getParam(token_param::CHR)){
                    m_chr=token.getParam(token_param::CHR);
                    writeHeader=true;
                }
                if (m_span!=ourspan){
                    m_span=ourspan;
                    writeHeader=true;
                }
            if (writeHeader)
            {
                *m_pOstream <<step<<m_chr<<"\t" <<span<<m_span<<'\n';
            }

                 *m_pOstream << ourstart << "\t" <<density << '\n';
           }
        else
            throw uWriter_missing_mandatory_param() << string_error("Token lacking mandatory param to write to Wig format. \n");
	}
	catch(param_not_found& e) {
		throw e;
	}
}

//DerivedRegister<uWriterWig> uWriterWig::reg("WIG");
} // End of namespace NGS
