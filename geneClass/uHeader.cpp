#include "uHeader.h"
namespace NGS {
/** \brief If the param exist, it returns it's value in string format. Otherwise returns param_not_found error.
 * \param header_param name the name of the param from header_param strongly type enum types.
 */
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

/** \brief To copy a uHeader object.
 */
uHeader& uHeader::operator=(uHeader const& assign_from) {
	if (this == &assign_from) return *this;
	m_params = assign_from.m_params;
	m_unformatedHeader = assign_from.m_unformatedHeader;
	return *this;
}

/** \brief Set the value of a param
 * \param const header_param& name: the name of the param. Must be a header_param type.
 * \param const std::string& value: the value of the param in string format. 
 */
void uHeader::_setParam(const header_param& name, const std::string& value){
	// TODO: Valiation and post-process of param?
	m_params[name] = value;
}

/** \brief Utility function to convert an header_param type to it's string counterpart.
 * \param const header_param& header: the header_param type to convert to string.
 */
std::string uHeader::_convertHeaderParamToString(const header_param& header) const {
	// TODO: Code operator << when we have some header_param
//	std::stringstream ss;
//	ss << token;
//	return ss.str();
	return "";
}
} // End of namespace NGS
