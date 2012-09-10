#include "uHeader.h"

std::string uHeader::getParam(header_param name) const{
	if(isParamSet(name)) {
		return (m_params.find(name)->second);
	}
	else {
		param_not_found e;
		std::string error = "Tried to getParam that is not set: " + _convertHeaderParamToString(name) + "\n";
		e << string_error(error);
		e << string_error(_convertHeaderParamToString(name));
		throw e;
	}
}

uHeader& uHeader::operator=(uHeader const& assign_from) {
	if (this == &assign_from) return *this;
	m_params = assign_from.m_params;
	m_unformatedHeader = assign_from.m_unformatedHeader;
	return *this;
}

void uHeader::_setParam(const header_param& name, const std::string& value){
	// TODO: Valiation and post-process of param?
	m_params[name] = value;
}

std::string uHeader::_convertHeaderParamToString(const header_param& header) const {
	// TODO: Code operator << when we have some header_param
//	std::stringstream ss;
//	ss << token;
//	return ss.str();
	return "";
}
