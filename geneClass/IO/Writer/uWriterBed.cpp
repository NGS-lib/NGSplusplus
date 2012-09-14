#include "uWriterBed.h"

namespace NGS {

/** \brief Destructor of uWriterBed
  */
uWriterBed::~uWriterBed() {
	if (m_dynamicStream == true) {
		delete m_pOstream;
	}
	m_pOstream = nullptr;
}

/** \brief Initialise the Bed writer with a file name.
  * \param const std::string& filename: name of the file to write output. 
  */
void uWriterBed::init(const std::string& filename) {
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


/** \brief Initialise the Bed writer with a stream
  * \param std::ostream* os: the stream to save the data to.
  */
void uWriterBed::init(std::ostream* os) {
	m_pOstream = os;
}

/** \brief Print the values of a token in Bed format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterBed::writeToken(const uToken& token) {
	try {
		std::string chr = token.getParam(token_param::CHR);	
		std::string start_pos = token.getParam(token_param::START_POS);	
		std::string end_pos = token.getParam(token_param::END_POS);	
		std::string seq_name = token.getParam(token_param::SEQ_NAME);	
		*m_pOstream << chr << '\t' << start_pos << '\t' << end_pos << '\t' << seq_name;
		/**< If score is set, we expect strand to be setted too */
		if (token.isParamSet(token_param::SCORE)) {
			std::string score = token.getParam(token_param::SCORE);	
			std::string strand = token.getParam(token_param::STRAND);	
			*m_pOstream << '\t' << score << '\t' << strand;
		}
		*m_pOstream << std::endl;
	}
	catch(param_not_found& e) {
		throw e;
	}
}

DerivedRegister<uWriterBed> uWriterBed::reg("BED");
} // End of namespace NGS
