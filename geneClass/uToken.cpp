#include "uToken.h"

//enum class token_param { CHR, START_POS, END_POS, STRAND, MAP_SCORE, PHRED_SCORE, CIGAR, SEQUENCE, SEQ_NAME, FLAGS };

/* \brief uToken constructor.
 * \param istream& paramList: a stream containing all the parameters in the format: <token_param>\t<value>\n<token_param>\t<value>\n...
 */
uToken::uToken(std::istream& paramList){
	char line[1024];	
	while (paramList.getline(line, 1024)) {
		std::stringstream ss;
		ss << line;
		/**< Fetch name */
		token_param name;
		try {
			ss >> name;
		}
		catch (invalid_token_param_error& e) {
			throw e;
		}
		/**< Fetch value */
		std::string value;
		ss >> value;
		try {
			_setParam(name, value);
		}
		catch (invalid_value_error& e) {
			throw e;
		}
	}
	/**< Check if uToken is in a valid state  */
	try {
		_validateToken();
	}
	catch (invalid_uToken_error& e) {
		throw e;
	}
}

//std::string uToken::getParam(token_param name) {
//	return m_params[name];
//}

/* \brief Set a param only if it's format is valid
 * \param token_param& name: type of param
 * \param const std::string& value: the value of the parameter
 */
void uToken::_setParam(token_param& name, std::string& value) {
	// TODO: What should we do in case the name is already setted? Overwrite silently, warning or error?
	if (_validateParam(name, value) == false) {
		invalid_value_throw e;
		e << invalid_value_error(value);
		throw e;
	}
	m_params[name] = value;
//	m_params.insert(make_pair(name, value));
}

/* \brief Validate the parameter. Some tests are general, other are specific for each parameter.
 * \param token_param& name: type of param
 * \param const std::string& value: the value of the parameter
 */
bool uToken::_validateParam(token_param& name, const std::string& value) const {
	bool isValid = true;
	/**< In every case, an empty string is an invalid value */
	if (value.size() == 0) {
		return false;
	}
	/**< Then we do basic check for every type of parameters */
	switch(name) {
//	case token_param::CHR:
//		isValid = _chrIsValid(value);
//		break;
	case token_param::START_POS:
	case token_param::END_POS:
		isValid = _posIsValid(value);
		break;
	case token_param::STRAND:
		isValid = _strandIsValid(value);
		break;
	case token_param::MAP_SCORE:
		isValid = _mapScoreIsValid(value);
		break;
	case token_param::PHRED_SCORE:
		isValid = _phredScoreIsValid(value);
		break;
	case token_param::CIGAR:
		isValid = _cigarIsValid(value);
		break;
	case token_param::SEQUENCE:
		isValid = _sequenceIsValid(value);
		break;
//	case token_param::SEQ_NAME:
//		isValid = _seqNameIsValid(value);
//		break;
	case token_param::FLAGS:
		isValid = _seqFlagsIsValid(value);
		break;
	defaut: break;
	}
	return isValid;
}

bool uToken::_validateToken() {
	/**< Validate start/end positions  */
	std::string str_start = getParam(token_param::START_POS);
	std::string str_end = getParam(token_param::END_POS);
	int start = atoi(str_start.c_str());
	int end = atoi(str_end.c_str());
	if (start > end) {
		std::string error = "Invalid START_POS/END_POS values. \n";
		error += "START_POS:" + str_start + ". END_POS: " + str_end + ".";
		_throwInvalidToken(error);
	}
	/**< Validate sequence/phred  */
	std::string str_phred = getParam(token_param::PHRED_SCORE);
	if (str_phred.size() != 0) {
		std::string str_sequence = getParam(token_param::SEQUENCE);
		if (str_sequence.size() != str_phred.size()) {
			std::stringstream seq_size;
			seq_size << str_sequence.size();
			std::stringstream phred_size;
			phred_size << str_phred.size();
			std::string error = "Sequence and phred score length does not match.\n";
			error += "Sequence size: " + seq_size.str() + ". Phred size: " + phred_size.str() + ".";
			_throwInvalidToken(error);
		}
	}
	/**< Validate sequence/cigar - Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ. (From SAM specifications) */
	std::string str_cigar = getParam(token_param::CIGAR);
	if (str_cigar.size() != 0) {
		std::string str_sequence = getParam(token_param::SEQUENCE);
		/**< Fetch all numerical values and sum them  */	
		size_t cigar_size = 0;
		std::stringstream ss;
		for (size_t i = 0; i < str_cigar.size(); i++) {
			if (_isDigit(str_cigar[i]) == true) {
				ss << str_cigar[i];
			}
			else {
				int value;
				switch(str_cigar[i]) {
				case 'M': 
				case 'I':
				case 'S':
				case '=':
				case 'X': ss >> value; break;
				default: value = 0; break;
				}
				cigar_size += value;
				ss.clear();
			}
		}
		/**< Check if sequence length and cigar sum match  */	
		if (str_sequence.size() != cigar_size) {
			std::stringstream seq_size_stream;
			seq_size_stream << str_sequence.size();
			std::stringstream cigar_size_stream;
			cigar_size_stream << cigar_size;
			std::string error = "Sequence length and sum of cigar M/I/S/=/X values do not match.\n";
			error += "Sequence length: " + seq_size_stream.str();
			error += ". Sum of cigar values: " + cigar_size_stream.str() + ".";
			_throwInvalidToken(error);
		}
	}
}

void uToken::_throwInvalidToken(const std::string& baseErrorMessage) const {
	invalid_uToken_throw e;
	e << invalid_uToken_error("Invalid token. " + baseErrorMessage);
	throw e;
}

// TODO: So far, it does not seem important to validate the CHR format, remove the commented section after tests are completed.
/* \brief Test if chr is valid.
 * \param const std::string& value: the value of the chr entry.
 */
//bool uToken::_chrIsValid(const std::string& value) const {
//	/**< Check if value starts with "chr" */
//	if (value[0] != 'c' && value[1] != 'h' && value[2] != 'r') {
//		return false;
//	}
//	/**< Check if next value is a positive int value */
//	std::stringstream ss;
//	ss << value.substr(3);
//	int n = 0;
//	ss >> n;
//	if (n <= 0) {
//		return false;
//	}
//	/**< Finally, check if there is garbage at the end of value */
//	return _isStreamEmpty(ss);
//}
// TODO: Old version, remove after testing is done
//	std::stringstream ss;
//	ss << value;
//	/**< Check if value starts with "chr" */
//	std::string chr;
//	ss >> chr;
//	if (chr != "chr") {
//		return false;
//	}
//	/**< Check if next value is a positive int value */
//	int n = 0;
//	ss >> n;
//	if (n <= 0) {
//		return false;
//	}
//	/**< Finally, check if there is garbage at the end of value */
//	return _isStreamEmpty(ss);
//}

/* \brief check if START_POS or END_POS is valid. 
 * \param const std::string& value: the value of the pos entry.
 */
bool uToken::_posIsValid(const std::string& value) const {
	std::stringstream ss;
	ss << value;
	/**< Check if the value is a long int greater than 0 */
	long int pos = 0;
	ss >> pos;
	if (pos <= 0) {
		return false;
	}
	/**< Finally, check if there is garbage at the end of value */
	return _isStreamEmpty(ss);
}

/* \brief Check if strand value is valid.
 * \param const std::string& value: the value of the strand entry.
 */
bool uToken::_strandIsValid(const std::string& value) const {
	/**< Check if size is equal to 1 */
	if (value.size() != 1) {
		return false;
	}
	/**< Check if value is either '+' or '-' */
	if (value != "+" || value != "-") {
		return false;
	}
	return true;
}
// TODO: Doc below this point may be copypasta, check to make sure it was updated correctly
/* \brief Check if map score value is valid (between 0 and 255 incl.).
 * \param const std::string& value: the value of the map score entry.
 */
bool uToken::_mapScoreIsValid(const std::string& value) const {
	std::stringstream ss;
	ss << value;
	/**< Check if value is between */
	int map = -1;
	ss >> map;
	if (map < 0 || map > 255) {
		return false;
	}
	/**< Finally, check if there is garbage at the end of value */
	return _isStreamEmpty(ss);
}

/* \brief Check phred score value is valid.
 * \param const std::string& value: the value of the phred score entry.
 */
bool uToken::_phredScoreIsValid(const std::string& value) const {
	/**< Check if value is a list of values (between 0 and 255 incl.). */
	for (size_t i = 0; i < value.size(); i++) {
//		const char* p = value[i].c_str();
		char p = value[i];
		if (p < 0 || p > 255) {
			return false;
		}
	}
}

/* \brief Check is sequence is valid.
 * \param const std::string& value: the value of the sequence entry.
 */
bool uToken::_sequenceIsValid(const std::string& value) const {
	/**< Check if the first value is '*' and the size is 1 */
	if (value[0]== '*' && value.size() > 1) {
		return false;
	}
	/**< Otherwise, check if all the characters are in a valid format */
	for (size_t i = 0; i < value.size(); i++) {
		switch (value[i]) {
		case 'a': break;
		case 'A': break;
		case 'c': break;
		case 'C': break;
		case 'g': break;
		case 'G': break;
		case 't': break;
		case 'T': break;
		case 'n': break;
		case 'N': break;
		case '.': break;
		default: return false;
		}
	}
}

// TODO: I don't think there is really anything worth checking for a sequence other once we know it's not empty
/* \brief Check if sequence name value is valid.
 * \param const std::string& value: the value of the sequence name entry.
 */
//bool uToken::_seqNameIsValid(const std::string& value) const {
//	std::stringstream ss;
//	ss << value;
//	/**< Check if value is either '+' or '-' */
//	/**< Finally, check if there is garbage at the end of value */
//	return _isStreamEmpty(ss);
//}

/* \brief Check if flag value is valid (between 0 and 65235 incl.).
 * \param const std::string& value: the value of the flag entry.
 */
bool uToken::_seqFlagsIsValid(const std::string& value) const {
	std::stringstream ss;
	ss << value;
	/**< Check if value is between 0 and 65235 incl. */
	int flags = -1;
	ss >> flags;
	if (flags < 0 || flags > 65235) {
		return false;
	}
	/**< Finally, check if there is garbage at the end of value */
	return _isStreamEmpty(ss);
}

/* \brief Check if cigar value is valid.
 * \param const std::string& value: the value of the cigar entry.
 */
bool uToken::_cigarIsValid(const std::string& value) const {
	/**< Check if the first value is '*' and the size is 1 */
	if (value[0] == '*' && value.size() > 1) {
		return false;
	}
	/**< First item of value has to be a digit.  */
	if (_isDigit(value[0]) != true) {
		return false;
	}
	/**< Make sure there is not two alphabetical character in a row */
	bool lastValue = true;
	for (size_t i = 0; i < value.size(); i++) {
		bool currentValue = _isDigit(value[i]);
		if (lastValue == false && currentValue == false) {
			return false;
		}
		lastValue = currentValue;
	}
	/**< Must end with alphabetical value */
	if (_isDigit(value[value.size()-1]) == true) {
		return false;
	}
	return true;
}

/* \brief Check if a value of a cigar score is a legal character.
 * \param char value: The value to test.
 */
bool uToken::_cigarValueIsValid(char value) const {
	switch(value) {
	case 'M': return true;
	case 'I': return true;
	case 'D': return true;
	case 'N': return true;
	case 'S': return true;
	case 'H': return true;
	case 'P': return true;
	case 'X': return true;
	case '=': return true;
	default: return false;
	}
}

/* \brief Check if a istream is empty.
 * \param const std::istream& stream: the istream to check.
 */
bool uToken::_isStreamEmpty(const std::istream& stream) const {
	if (stream.rdbuf()->in_avail() != 0) {
		return false;
	}
	return true;
}
