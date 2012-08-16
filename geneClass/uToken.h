#ifndef UTOKEN_H_INCLUDED
#define UTOKEN_H_INCLUDED

#include <map>
#include <string>
#include <stdexcept>
#include "boost/exception/all.hpp"

/**< List of param is hard coded as strongly typed enum for extra safety. */
/**< This list has to be updated for every new param */
enum class token_param { CHR, START_POS, END_POS, STRAND, MAP_SCORE, PHRED_SCORE, CIGAR, SEQUENCE, SEQ_NAME, FLAGS };

/**< uToken class, to bridge parser and the library's class */
/**< This is the class that takes care of data validation */
/**< All the data is saved in a map in string format */
class uToken {
public:
	uToken() { };
	void setParam(token_param::token_param& name, const std::string& value);
	std::string getParam(token_param::token_param& name) const;

private:
	std::map<token_param, std::string> m_params;
	void _validateParam(token_param::token_param& name, const std::string& value);

} // End of class Token


// uToken exceptions
struct uToken_exception_base : virtual std::exception, virtual boost::exception {};
struct value_throw : virtual uToken_exception_base{};
typedef boost::error_info<struct invalid_param_info, std::string> invalid_value_error;

#endif // UTOKEN_H_INCLUDED
