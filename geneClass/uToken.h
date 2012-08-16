#ifndef UTOKEN_H_INCLUDED
#define UTOKEN_H_INCLUDED

#include <map>
#include <string>
#include <sstream>
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
	uToken(istream& paramList);
	std::string getParam(token_param::token_param& name) const;

private:
	std::map<token_param, std::string> m_params;
	void _setParam(token_param::token_param& name, const std::string& value);
	void _validateParam(token_param::token_param& name, const std::string& value);
	bool _chrIsValid(const std::string& name) const;
	bool _posIsValid(const std::string& name) const;
	bool _strandIsValid(const std::string& name) const;

} // End of class Token

/**< Overloading of stream operator for token_param */
inline std::ostream& operator<<(std::ostream & Str, token_param::token_param name) {
	switch (name) {
	case CHR: return Str << "CHR";
	case START_POS: return Str << "START_POS";
	case END_POS: return Str << "END_POS";
	case STRAND: return Str << "STRAND";
	case MAP_SCORE: return Str << "MAP_SCORE";
	case PHRED_SCORE: return Str << "PHRED_SCORE";
	case CIGAR: return Str << "CIGAR";
	case SEQUENCE: return Str << "SEQUENCE";
	case SEQ_NAME: return Str << "SEQ_NAME";
	case FLAGS: return Str << "FLAGS";
	default: return Str << (int) name;
	}
}

inline istream& operator >>(istream &is, token_param::token_param& name) {
	is >> name;
	return is;
}

/**<  uToken exceptions */
struct uToken_exception_base : virtual std::exception, virtual boost::exception {};
struct invalid_value_throw : virtual uToken_exception_base{};
typedef boost::error_info<struct invalid_param_info, std::string> invalid_value_error;

#endif // UTOKEN_H_INCLUDED
