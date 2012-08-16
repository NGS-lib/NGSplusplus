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
		/**< Fetch value */
	}

}

/* \brief Fetch a param. Throw param_not_found if the param does not exist.
 * \param token_param& name: the name of the param we wish to get.
 */
std::string uToken::getParam(token_param& name) const {

}

/* \brief Set a param only if it's format is valid
 * \param token_param& name: type of param
 * \param const std::string& value: the value of the parameter
 */
void uToken::_setParam(token_param& name, std::string& value) {
	// TODO: What should we do in case the name is already setted? Overwrite silently, warning or error?
	try {
		_validateParam(name, value);
	}
	catch (invalid_value_error& e) {
		throw e;		
	}
	m_params.insert(make_pair(name, value));
}

/* \brief Validate the parameter. Some tests are general, other are specific for each parameter.
 * \param token_param& name: type of param
 * \param const std::string& value: the value of the parameter
 */
void uToken::_validateParam(token_param& name, const std::string& value) {
	/**< In every case, an empty string is an invalid value */
	if (value.size() == 0) {
		invalid_value_throw e;
		e << invalid_value_error(value);
		throw e;
	}
	/**< Then we do basic check for every type of parameters */
	switch(name) {
	case token_param::CHR:
		break;
	case token_param::START_POS:
		break;
	case token_param::END_POS:
		break;
	case token_param::STRAND:
		break;
	case token_param::MAP_SCORE:
		break;
	case token_param::PHRED_SCORE:
		break;
	case token_param::CIGAR:
		break;
	case token_param::SEQUENCE:
		break;
	case token_param::SEQ_NAME:
		break;
	case token_param::FLAGS:
		break;
	}
}
