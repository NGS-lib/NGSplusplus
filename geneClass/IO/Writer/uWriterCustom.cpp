#include "uWriterCustom.h"

namespace NGS {

/** \brief Destructor of uWriterCustom
  */
uWriterCustom::~uWriterCustom() {
	if (m_dynamicStream == true) {
		delete m_pOstream;
	}
	m_pOstream = nullptr;
}

/** \brief Initialise the Custom writer with a file name.
  * \param const std::string& filename: name of the file to write output. 
  */
void uWriterCustom::init(const std::string& filename) {
	std::filebuf fb;
	fb.open (filename.c_str(),std::ios::out);
	if (!fb.is_open()) {  
		std::string error = "Error opening file: " + filename;
		throw std::runtime_error(error.c_str());
	}
	else {  
		std::ostream* os = new std::ostream(&fb);
		m_pOstream = os;
	}
	m_dynamicStream = true;
}

/** \brief Initialise the Custom writer with a stream
  * \param std::ostream* os: the stream to save the data to.
  */
//void uWriterCustom::init(std::ostream* os, const std::vector<std::string>& fieldsNames) {
void uWriterCustom::init(std::ostream* os) {
	if (os != nullptr && os->good() == true) {
		m_pOstream = os;
	}
	else {
		throw std::runtime_error("Invalid stream.");
	}
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
			*m_pOstream << '\t';
		}
	}
	*m_pOstream << std::endl;
}

DerivedRegister<uWriterCustom> uWriterCustom::reg("CUSTOM");
} // End of namespace NGS
