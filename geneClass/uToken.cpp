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
//	if (_validateToken() == false) {
//
//	}
}

/* \brief Fetch a param. Throw param_not_found if the param does not exist.
 * \param token_param& name: the name of the param we wish to get.
 */
std::string uToken::getParam(token_param& name) const {
	return "asdf";
}

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
	m_params.insert(make_pair(name, value));
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
	/**< Check if length of value is an even number */
	if (value.size() % 2  != 0) {
		return false;
	}
	/**< Check if the first value is '*' and the size is 1 */
	if (value[0] == '*' && value.size() > 1) {
		return false;
	}
	/**< Check if even positions in the value are numbers */
	for (size_t i = 0; i < value.size(); i = i + 2) {
		if (value[i] < '0' || value[i] > '9') {
			return false;
		}
	}
	/**< Check if odd positions in the value are valid characters */
	for (size_t i = 1; i < value.size(); i = i + 2) {
		switch(value[i]) {
		case 'M': break;
		case 'I': break;
		case 'D': break;
		case 'N': break;
		case 'S': break;
		case 'H': break;
		case 'P': break;
		case 'X': break;
		default: return false;
		}
	}
	return true;
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
