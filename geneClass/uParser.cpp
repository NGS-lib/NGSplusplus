#include "uParser.h"
#include "uGeneException.h"
/** \brief uParser constructor with filename
 * \param std::string filename: Name of the file to load
 * \param file_type type: Currently supported formats: BED, SAM
 */
uParser::uParser(const std::string& filename, file_type type, bool header) {
	std::ifstream* ifs = new std::ifstream(filename.c_str(), std::ifstream::in);
	if (!ifs->is_open()) {
		std::string error = "Error opening file: " + filename;
		throw std::runtime_error(error.c_str());
	}
	else {
		m_pIstream = ifs;
	}
	m_fileType = type;
	m_header = header;
	if (m_header == false) {
		m_firstToken = false;
	}
	m_dynamicStream = true;
}

/** \brief uParser constructor with istream
 * \param std::istream* stream: stream to load the data from
 * \param file_type type: Currently supported formats: BED, SAM
 */
uParser::uParser(std::istream* stream, file_type type, bool header) {
	m_pIstream = stream;
	m_fileType = type;
	m_header = header;
	if (m_header == false) {
		m_firstToken = false;
	}
	m_dynamicStream = false;
}

/** \brief uParser constructor with custom file type
 * \param const std::vector<string> columnNames: The name of every field in the custom format, in the SAME ORDER as they appear in the file. Must be of the string version of token_param (see uToken.h). Mandatory fields are: CHR, START_POS and END_POS.
 * \param char delimiter: The delimiter between each field.
 */
uParser::uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header, char delimiter) {
	/**< Check if filename is valid, then open it */
	std::ifstream* ifs = new std::ifstream(filename.c_str(), std::ifstream::in);
	if (!ifs->is_open()) {
		std::string error = "Error opening file: " + filename;
		throw std::runtime_error(error.c_str());
	}
	else {
		m_pIstream = ifs;
	}
	/**< Check if fields are in a valid format */
	try {
		_customParserValidateFields(fieldsNames);
		_customParserCopyFields(fieldsNames);
	}
	catch(customParser_missing_mandatory_values& e) {
		throw e;
	}
	/**< Set other parameters */
	m_fileType = file_type::CUSTOM;
	m_delimiter = delimiter;
	m_header = header;
	if (m_header == false) {
		m_firstToken = false;
	}
	m_dynamicStream = true;
}

uParser::uParser(std::istream* stream, const std::vector<std::string>& fieldsNames, bool header, char delimiter) {
	/**< Check if fields are in a valid format */
	try {
		_customParserValidateFields(fieldsNames);
		_customParserCopyFields(fieldsNames);
	}
	catch(customParser_missing_mandatory_values& e){
		throw e;
	}
	/**< Set other parameters */
	m_pIstream = stream;
	m_fileType = file_type::CUSTOM;
	m_delimiter = delimiter;
	m_header = header;
	if (m_header == false) {
		m_firstToken = false;
	}
	m_dynamicStream = false;
}

/** \brief uParser destructor
 */
uParser::~uParser() {
	if (m_dynamicStream == true) {
		delete m_pIstream;
	}
	m_pIstream = NULL;
}

/** \brief Generic function to fetch an entry from file, will adapt itself to file_type given in constructor
 */
uToken uParser::getNextEntry() {
	try {
		switch(m_fileType) {
		case file_type::BED:
			try {
				uToken token = _getNextEntryBed(); return token;
			}
			catch(invalid_uToken_throw& e) {
				#ifdef DEBUG
				std::string trace = fetchStringError(e);
				std::cerr << "Invalid uToken: " << trace << std::endl;
				#endif
				throw e;
			}
			catch(end_of_file_throw& e) {
				throw e;
			}
			break;

		case file_type::SAM:
			try {
				uToken token = _getNextEntrySam(); return token;
			}
			catch(invalid_uToken_throw& e) {
				#ifdef DEBUG
				std::string trace = fetchStringError(e);
				std::cerr << "Invalid uToken: " << trace << std::endl;
				#endif
				throw e;
			}
			break;

		case file_type::CUSTOM:
			try {
				uToken token = _getNextEntryCustom(); return token;
			}
			catch(invalid_uToken_throw& e) {
				#ifdef DEBUG
				std::string trace = fetchStringError(e);
				std::cerr << "Invalid uToken: " << trace << std::endl;
				#endif
				throw e;
			}

		default:
			throw uParser_exception_base()<< string_error("Invalid fileType in getNextEntry case");
			break;
		}
	}
	catch(invalid_uToken_throw &e){
		throw e;
	}
}

// TODO: class uHeader to return from this function. Derived from uToken?
/** \brief Fetch header differently based on file type
 */
void uParser::_fetchHeader() {
	/**< When header is set as true for BED and CUSTOM, we do nothing here. */
	/**< Instead, the code won't throw error for bad token until first valid token. */
	switch(m_fileType) {
	case file_type::BED:
	case file_type::CUSTOM: break;
	// TODO: Write header parser for SAM format
	case file_type::SAM: break;
	default: break;
	}
}

/** \brief Specific loader for BED file (See genome.ucsc.edu/FAQ/FAQformat.html#format1 for bed description)
 * \return uToken: If all the parameters in the entry are valid a uToken object is returned.
 */
uToken uParser::_getNextEntryBed() {
	do {
		char line[4096];
		if (m_pIstream->getline(line, 4096)) {
			/**< We start by fetching the infos from the line */
			std::stringstream ss;
			ss << line;
			std::string chr;
			std::string start_pos;
			std::string end_pos;
			std::string score;
			std::string seq_name;
			std::string strand;
			std::stringstream token_infos;
			ss >> chr >> start_pos >> end_pos >> seq_name >> score >> strand;
			token_infos << "CHR\t" << chr << "\n";
			token_infos << "START_POS\t" << start_pos << "\n";
			token_infos << "END_POS\t" << end_pos << "\n";
			token_infos << "SEQ_NAME\t" << seq_name << "\n";

			/**< If there was no strand info, we don't add an empty string */
			if (strand.size() != 0) {
				token_infos << "STRAND\t" << strand << "\n";
			}
			/**< We try to create a token with the infos that were fetched from the line */
			/**< If it doesn't work andit's the first token, we don't throw an error. Instead, we try again with another line */
			try {
				uToken token(token_infos);
				m_firstToken = false;
				return token;
			}
			catch(invalid_uToken_throw& e) {
				if (m_firstToken != true) {
					throw e;
				}
			}
			catch(invalid_value_throw& e) {
				if (m_firstToken != true) {
					throw e;
				}
			}
		}
		else {
			#ifdef DEBUG
			std::cerr << "Reached end of file." << std::endl;
			#endif
			end_of_file_throw e;
			e << string_error("Reached end of file.");
			throw e;
		}
	}while (m_firstToken == true);

std::cerr <<"Fatal error in _getNextEntryCustom(), should not reach here" <<std::endl;
abort();
}

/** \brief Specific loader for SAM file (See samtools.sourceforge.net for SAM description)
 * \return uToken: If all the parameters in the entry are valid a uToken object is returned.
 */
uToken uParser::_getNextEntrySam() {
	try {
		char line[4096];
		if (m_pIstream->getline(line, 4096)) {
			std::stringstream ss;
			ss << line;

			//String to test for int
			std::string flag;
			std::string start_pos;
			std::string MAPQual;
			std::string pNext;
			std::string Tlen;

			std::string chr;

			std::string end_pos;
			std::string score;
			std::string seq_name;
			std::string qual;
			std::string seq;
			std::string cigar;
			std::string RNext;
			std::stringstream token_infos;

			ss >> seq_name >> flag >> chr >> start_pos >> MAPQual >> cigar>>RNext>>pNext>>Tlen>>seq>>qual;

			token_infos << "CHR\t" << chr << "\n";
			token_infos << "START_POS\t" << start_pos << "\n";
			//token_infos << "END_POS\t" << end_pos << "\n";
			token_infos << "FLAG\t" << flag << "\n";
			token_infos << "SEQ_NAME\t" << seq_name << "\n";
			token_infos << "MAP_SCORE\t" << MAPQual << "\n";
			token_infos << "SEQUENCE\t" << seq << "\n";
			token_infos << "CIGAR\t" << cigar << "\n";
			token_infos << "PHRED_SCORE\t" << qual << "\n";

			uToken token(token_infos);
			return token;
		}
		else {
			#ifdef DEBUG
			std::cerr << "Reached end of file." << std::endl;
			#endif
			end_of_file_throw e;
			e << string_error("Reached end of file.");
			throw e;
		}
	}
	catch(invalid_uToken_throw& e) {
		throw e;
	}
}

uToken uParser::_getNextEntryCustom() {
	do {
		char line[4096];
		if (m_pIstream->getline(line, 4096)) {
			/**< We start by fetching the infos in the line */
			std::stringstream token_infos;
			char* current;
			current = strtok(line, &m_delimiter);
			for(size_t i = 0; i < m_customFieldNames.size(); i++) {
				if (m_customFieldNames[i] != "NA") {
					token_infos << m_customFieldNames[i] << "\t" << current << "\n";
				}
				current = strtok(NULL, &m_delimiter);
			}
			/**< Then we try to create a token with that info */
			/**< If it doesn't work andit's the first token, we don't throw an error. Instead, we try again with another line */
			try {
				uToken token(token_infos);
				m_firstToken = false;
				return token;
			}
			catch(invalid_uToken_throw& e) {
				if (m_firstToken != true) {
					throw e;
				}
			}
			catch(invalid_value_throw& e) {
				if (m_firstToken != true) {
					throw e;
				}
			}
		}
		else {
			#ifdef DEBUG
			std::cerr << "Reached end of file." << std::endl;
			#endif
			end_of_file_throw e;
			e << string_error("Reached end of file.");
			throw e;
		}
	} while(m_firstToken == true);

std::cerr <<"Fatal error in _getNextEntryCustom(), should not reach here" <<std::endl;
abort();
}

void uParser::_customParserValidateFields(const std::vector<std::string>& fieldsNames) {
	/**< Must have at least 2 fields */
	if (fieldsNames.size() < 3) {
		customParser_missing_mandatory_values e;
		e << string_error("Custom file fields description is too short, must have at least 3 values: CHR and START_POS and a way to infer END_POS.\n");
		throw e;
	}

	/**< Check if mandatory fields are present */
	if (!_paramExists("CHR", fieldsNames)) {
		customParser_missing_mandatory_values e;
		e << string_error("Mandatory field is missing: CHR\n");
		throw e;
	}
	if (!_paramExists("START_POS", fieldsNames)) {
		customParser_missing_mandatory_values e;
		e << string_error("Mandatory field is missing: START_POS\n");
		throw e;
	}

	/**< We need to be able to infer END_POS either directly or indirectly (with a sequence or cigar score) */
	if (!_paramExists("END_POS", fieldsNames) && !_paramExists("SEQUENCE", fieldsNames) && !_paramExists("CIGAR", fieldsNames)) {
		customParser_missing_mandatory_values e;
		e << string_error("We must be able to infer END_POS directly or indirectly (with SEQUENCE or CIGAR)\n");
		throw e;
	}
}

bool uParser::_paramExists(const std::string& name, const std::vector<std::string>& list) const {
	std::vector<std::string>::const_iterator it;
	it = find(list.begin(), list.end(), name);
	if (it == list.end()) {
		return false;
	}
	return true;
}

void uParser::_customParserCopyFields(const std::vector<std::string>& fieldsNames) {
	/**< Only copy fields that have a matching token_param value, otherwise add NA */
	for(size_t i = 0; i < fieldsNames.size(); i++) {
		if (uToken::checkParam(fieldsNames[i]) == true) {
			m_customFieldNames.push_back(fieldsNames[i]);
		}
		else {
			m_customFieldNames.push_back("NA");
		}
	}
}
