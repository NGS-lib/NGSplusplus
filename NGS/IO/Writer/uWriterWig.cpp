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
                long long int ourstart= std::stoll(token.getParam(token_param::START_POS));
                long long int ourend= std::stoll(token.getParam(token_param::END_POS));
                int ourspan = ( ourend-ourstart+1);
                float density =1;
                if (token.isParamSet(m_scoreType))
                {
                 /**< Should not be using catch as a control flow method, but will to for now */
                 /**< Establish standard for bed conversion */
                 try {
                    density=std::stof(token.getParam(m_scoreType));
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
                *m_pOstream <<step<<m_chr<<"\t" <<span<<m_span<<endl;
            }

                 *m_pOstream << ourstart << "\t" <<density << endl;
           }
        else
            throw uWriter_missing_mandatory_param() << string_error("Token lacking mandatory param to write to Wig format. \n");

	}
	catch(param_not_found& e) {
		throw e;
	}
}

DerivedRegister<uWriterWig> uWriterWig::reg("WIG");
} // End of namespace NGS
