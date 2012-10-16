#include "uHeader.h"
namespace NGS {
/** \brief If the param exist, it returns it's value in string format. Otherwise returns param_not_found error.
 * \param header_param name the name of the param from header_param strongly type enum types.
 */

bool uHeader::_validateChrSize(std::string sizeString)
{
    try {
        int value = std::stoi(sizeString);
        /**< Chr size must be above 0 */
        if (value <= 0)
            return false;
        return true;
    }
    catch(...){
        return false;
    }
}

bool uHeader::_validateChrList(std::string chrString)
{
    /**< Make sur our name is unique */
    if (isParamSet(header_param::CHR_LIST))
    {
       auto listVector = getParamVector(header_param::CHR_LIST);
       auto it = find (listVector.begin(), listVector.end(), chrString);
       if (it!=listVector.end())
        return false;
    }

    return true;
}


 uHeader::uHeader(){

/**< Register our functions */
    /**< Validate functions */
    validate_func_map[header_param::CHR_LIST]= &uHeader::_noValidate;
    validate_func_map[header_param::CHR_SIZES]=&uHeader::_noValidate;
    /**< Post processing function */
    post_func_map[header_param::CHR_LIST]=&uHeader::_noPost;
    post_func_map[header_param::CHR_SIZES]=&uHeader::_noPost;

 }

std::vector<std::string> uHeader::getParamVector(header_param name) const{
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

std::string uHeader::getParam(header_param name) const{
	if(isParamSet(name)) {
		return (m_params.find(name)->second.at(0));
	}
	else {
		param_not_found e;
		std::string error = "Tried to getParam that is not set: " + _convertHeaderParamToString(name) + "\n";
		e << string_error(error);
		e << string_error(_convertHeaderParamToString(name));
		throw e;
	}
}

 //auto begin()->decltype(ExpMap.begin()){return ExpMap.begin();};
//auto uHeader::getParam(header_param name)->( std::tuple_element<0,  m_params.find(name)->second>)
//{
//    return std::get<0>(m_params.find(name)->second)
//}

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
	if (m_params[name].size()==0)
          m_params[name].resize(1);

    m_params[name].at(0) = value;
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

/** \brief Operate any necessary post_processing on our values
 *
 * \param name header_param& Type of parameter
 * \param value const std::string& Value of parameter
 * \return void
 *
 */
void uHeader::_postProcessParam(const header_param& name, const std::string& value) {
	post_func_map[name](this,value);
}


} // End of namespace NGS
