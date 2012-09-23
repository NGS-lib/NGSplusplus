#include "uWriterCustom.h"

namespace NGS {

uWriterCustom::uWriterCustom() {
	m_delimiter = '\t';
	m_fieldsNames.push_back("CHR");
	m_fieldsNames.push_back("START_POS");
	m_fieldsNames.push_back("END_POS");
	m_fieldsNames.push_back("SEQ_NAME");
}

/** \brief Print the values of a token in Custom format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterCustom::writeToken(const uToken& token) {
	std::vector<std::string> values;
	for (size_t i = 0; i < m_fieldsNames.size(); i++) {
		if (uToken::checkParam(m_fieldsNames[i]) == true) {
			std::stringstream ss;
			ss << m_fieldsNames[i];
			token_param param;
			try {
				ss >> param;
				if (token.isParamSet(param)) {
					values.push_back(token.getParam(param));
				}
				else {
					values.push_back(".");
				}
			}
			catch (invalid_token_param_throw& e) {
				values.push_back(".");
			}
		}
		else {
			values.push_back(".");
		}
		*m_pOstream << values[i];
		if (i != m_fieldsNames.size() - 1) {
			*m_pOstream << m_delimiter;
		}
	}
	*m_pOstream << std::endl;
}

DerivedRegister<uWriterCustom> uWriterCustom::reg("CUSTOM");

} // End of namespace NGS
