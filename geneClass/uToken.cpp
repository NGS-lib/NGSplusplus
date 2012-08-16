#include "uToken.h"

//enum class token_param { CHR, START_POS, END_POS, STRAND, MAP_SCORE, PHRED_SCORE, CIGAR, SEQUENCE, SEQ_NAME, FLAGS };

void uToken::setParam(token_param::token_param& name, const std::string& value) {
	// TODO: What should we do in case the name is already setted? Overwrite silently, warning or error?
	_validateParam(name, value);
	m_params(name) = value;
}

std::string uToken::getParam(token_param::token_param& name) const {

}

void uToken::_validateParam(token_param::token_param& name, const std::string& value) {
	
}
