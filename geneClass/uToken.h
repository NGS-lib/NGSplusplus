#ifndef UTOKEN_H_INCLUDED
#define UTOKEN_H_INCLUDED

#include <iostream>
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
	uToken(std::istream& paramList);
	std::string getParam(token_param& name) const;
//	friend std::ostream& operator<<(std::ostream& Str, token_param name);

private:
	std::map<token_param, std::string> m_params;
	void _setParam(token_param& name, std::string& value);
	void _validateParam(token_param& name, const std::string& value);
	bool _chrIsValid(const std::string& name) const { return true; };
	bool _posIsValid(const std::string& name) const { return true; };
	bool _strandIsValid(const std::string& name) const { return true; };
	bool _mapScoreIsValid(const std::string& name) const { return true; };
	bool _phredScoreIsValid(const std::string& name) const { return true; };

}; // End of class Token

/**< Overloading of stream operator for token_param */
inline std::ostream & operator<<(std::ostream& Str, token_param name) {
	switch (name) {
	case token_param::CHR: return Str << "CHR";
	case token_param::START_POS: return Str << "START_POS";
	case token_param::END_POS: return Str << "END_POS";
	case token_param::STRAND: return Str << "STRAND";
	case token_param::MAP_SCORE: return Str << "MAP_SCORE";
	case token_param::PHRED_SCORE: return Str << "PHRED_SCORE";
	case token_param::CIGAR: return Str << "CIGAR";
	case token_param::SEQUENCE: return Str << "SEQUENCE";
	case token_param::SEQ_NAME: return Str << "SEQ_NAME";
	case token_param::FLAGS: return Str << "FLAGS";
	default: return Str << (int) name;
	}
}

inline std::istream& operator>>(std::istream &is, token_param& name) {
	std::string token;
	is >> token;
	if (token == "CHR") name = token_param::CHR;
	else if (token == "START_POS") name = token_param::START_POS;
	else if (token == "END_POS") name = token_param::END_POS;
	else if (token == "STRAND") name = token_param::STRAND;
	else if (token == "MAP_SCORE") name = token_param::MAP_SCORE;
	else if (token == "PHRED_SCORE") name = token_param::PHRED_SCORE;
	else if (token == "CIGAR") name = token_param::CIGAR;
	else if (token == "SEQUENCE") name = token_param::SEQUENCE;
	else if (token == "SEQ_NAME") name = token_param::SEQ_NAME;
	else if (token == "FLAGS") name = token_param::FLAGS;
	return is;
}

/**<  uToken exceptions */
struct uToken_exception_base : virtual std::exception, virtual boost::exception {};
struct invalid_value_throw : virtual uToken_exception_base{};
typedef boost::error_info<struct invalid_param_info, std::string> invalid_value_error;

#endif // UTOKEN_H_INCLUDED
