#include "uParser.h"

/* \brief uParser constructor with filename
 * \param std::string filename: Name of the file to load
 * \param file_type type: Currently supported formats: BED
 */
uParser::uParser(const std::string& filename, file_type type) {
	std::ifstream* ifs = new std::ifstream(filename.c_str(), std::ifstream::in);
	if (!ifs->is_open()) {
		std::string error = "Error opening file: " + filename;
		throw std::runtime_error(error.c_str());
	}
	else {
		m_pIstream = ifs;
	}
	m_fileType = type;
	m_dynamicStream = true;
}

/* \brief uParser constructor with istream
 */
uParser::uParser(std::istream* stream, file_type type) {
	m_pIstream = stream;
	m_fileType = type;
	m_dynamicStream = false;
}

/* \brief uParser destructor
 */
uParser::~uParser() {
	if (m_dynamicStream == true) {
		delete m_pIstream;
	}
	m_pIstream = NULL;
}

/* \brief Generic function to fetch an entry from file, will adapt itself to file_type given in constructor
 */
uToken uParser::getNextEntry() {
	switch(m_fileType) {
	case file_type::BED: 
		try {
			uToken token = _getNextEntryBed(); return token;
		}
		catch(invalid_uToken_error& e) {
			// TODO: what should we do in case of error?
			// TODO: cerr warning?
			throw e;
		}
		catch(end_of_file_error& e) {
			// TODO: What do we do with eof? return ptr instead? smart ptr?
			throw e;
		}
		break;
	default: break;
	}
}

/* \brief Specific loader for BED file (See genome.ucsc.edu/FAQ/FAQformat.html#format1 for bed description)
 * \return uToken: If all the parameters in the entry are valid a uToken object is returned.
 */
uToken uParser::_getNextEntryBed() {
	if (!m_pIstream->eof()) {
		char line[4096];
		m_pIstream->getline(line, 4096);
		std::stringstream ss;
		ss << line;
		std::string chr;
		std::string start_pos;
		std::string end_pos;
		std::string score;
		std::string strand;
		std::stringstream token_infos;
		ss >> chr >> start_pos >> end_pos >> score >> strand;
		token_infos << "CHR\t" << chr << "\n";
		token_infos << "START_POS\t" << start_pos << "\n";
		token_infos << "END_POS\t" << end_pos;
		/**< If there was no strand info, we don't add an empty string */
		if (strand.size() != 0) {
			token_infos << "STRAND\t" << strand << "\n";
		}
		try {
			uToken token(token_infos);
		}
		catch(invalid_uToken_error& e) {
			// TODO: What to do with error?
			throw e;
		}
	}
	else {
		// TODO: What do we do when it is end of file?
		end_of_file_throw e;
		e << end_of_file_error("Reached end of file.");
		throw e;
	}
}

