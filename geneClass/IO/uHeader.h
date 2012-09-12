#ifndef UHEADER_H_INCLUDED
#define UHEADER_H_INCLUDED

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include "uGeneException.h"
#include "uGeneException.h"
namespace NGS {
/**< List of param is hard coded as strongly typed enum for extra safety. */
/**< This list has to be updated for every new param */
// TODO: Add param as they are needed. Remove TMP when there is one valid param, it's there to avoid compiling error.
enum class header_param { TMP };

/**< uHeader class, to keep track of information in header, formated or not */
// TODO: Data validation?
/**< This is the class that takes care of data validation */
/**< All the data is saved in a map in string format */
class uHeader {
public:
	/** \brief Default constructor.
	  */
	uHeader() {};
	void setUnformatedHeader(const std::string& unformatedHeader) { m_unformatedHeader = unformatedHeader; };
	/** \brief Fetch a param. Throw param_not_found if the param does not exist.
	 * \param header_param& name: the name of the param we wish to get.
	 */
	std::string getParam(header_param name) const;

	/** \brief Return the exact same header that was in the input data.
	  */
	std::string getUnformatedHeader() const { return m_unformatedHeader;};
	/** \brief Check if there is a value associated with a given param.
          * \param header_param& name: name of the param to check.
	  */
	bool isParamSet(const header_param& name) const { return m_params.count(name); };
	uHeader& operator=(uHeader const& assign_from);

	/** \brief Check if a string has a corresponding header_param value
	 * \param const std::string& param: The param to check.
	 * \return True if there is a matching header_param value, otherwise false.
	 */
	static inline bool checkParam(const std::string& param) { //TODO
		return true;
//		return (param == "CHR"
//		     || param == "START_POS"
//		     || param == "END_POS"
//		     || param == "STRAND"
//		     || param == "MAP_SCORE"
//		     || param == "PHRED_SCORE"
//		     || param == "CIGAR"
//		     || param == "SEQUENCE"
//		     || param == "SEQ_NAME"
//		     || param == "FLAGS");
	}

private:
	std::map<header_param, std::string> m_params={};
	std::string m_unformatedHeader="";
	void _setParam(const header_param& name, const std::string& value);
	std::string _convertHeaderParamToString(const header_param& header) const;

//	void _postProcessParam(const header_param& name, const std::string& value);
//	bool _validateParam(const header_param& name, const std::string& value) const;
//	void _validateHeader();
//	void _checkMandatoryValues() const;
//	void _throwInvalidHeader(const std::string& baseErrorMessage) const;

//	bool _isDigit(char value) const { return std::isdigit(value);}//(value >= '0' && value <= '9'); }
//	bool _isStreamEmpty(const std::istream& stream) const;

}; // End of class Header

/**< Overloading of stream operator for header_param */
//inline std::ostream & operator<<(std::ostream& Str, header_param name) {
//	switch (name) {
//	case header_param::CHR: return Str << "CHR";
//	case header_param::START_POS: return Str << "START_POS";
//	case header_param::END_POS: return Str << "END_POS";
//	case header_param::STRAND: return Str << "STRAND";
//	case header_param::MAP_SCORE: return Str << "MAP_SCORE";
//	case header_param::PHRED_SCORE: return Str << "PHRED_SCORE";
//	case header_param::CIGAR: return Str << "CIGAR";
//	case header_param::SEQUENCE: return Str << "SEQUENCE";
//	case header_param::SEQ_NAME: return Str << "SEQ_NAME";
//	case header_param::FLAGS: return Str << "FLAGS";
//	default: return Str << (int) name;
//	}
//}
//
//inline std::istream& operator>>(std::istream &is, header_param& name) {
//	std::string header;
//	is >> header;
//	if (header == "CHR") name = header_param::CHR;
//	else if (header == "START_POS") name = header_param::START_POS;
//	else if (header == "END_POS") name = header_param::END_POS;
//	else if (header == "STRAND") name = header_param::STRAND;
//	else if (header == "MAP_SCORE") name = header_param::MAP_SCORE;
//	else if (header == "PHRED_SCORE") name = header_param::PHRED_SCORE;
//	else if (header == "CIGAR") name = header_param::CIGAR;
//	else if (header == "SEQUENCE") name = header_param::SEQUENCE;
//	else if (header == "SEQ_NAME") name = header_param::SEQ_NAME;
//	else if (header == "FLAGS") name = header_param::FLAGS;
//	else {
//		invalid_header_param_throw e;
//		e << string_error(header);
//		throw e;
//	}
//	return is;
//}
} // End of namespace NGS
#endif // UHEADER_H_INCLUDED
