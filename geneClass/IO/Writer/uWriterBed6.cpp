#include "uWriterBed6.h"

namespace NGS {

/** \brief Print the values of a token in Bed format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterBed6::writeToken(const uToken& token) {
	try {
		uWriterBed::writeToken(token);
		/**< If score or strand are not set, we replace them with the value "." */
		std::string score;
		std::string strand;
		if (token.isParamSet(token_param::SCORE)) {
			score = token.getParam(token_param::SCORE);	
		}
		else {
			score = ".";
		}
		/**< If score is not set, but strand is, we replace score by "." */
		if (token.isParamSet(token_param::STRAND)) {
			strand = token.getParam(token_param::STRAND);
		}
		else {
			strand = ".";
		}
		*m_pOstream << '\t' << score << '\t' << strand;
	}
	catch(param_not_found& e) {
		throw e;
	}
	*m_pOstream << std::endl;
}

DerivedRegister<uWriterBed6> uWriterBed6::reg("BED6");
} // End of namespace NGS
