#ifndef UPARSER_H_INCLUDED
#define UPARSER_H_INCLUDED

#include <string>
#include <fstream>
#include "uToken.h"

enum class file_type{ BED, SAM };

class uParser {
public:
	uParser(std::string filename, file_type type = BED);
	uToken getNextEntry();

private:
	std::ifstream m_ifstream;
	token_param m_fileType;
	uToken _getNextEntryBed();
};

/**<  uParser exceptions */
struct uParser_exception_base : virtual std::exception, virtual boost::exception {};

struct invalid_entry_throw : virtual uToken_exception_base{};
typedef boost::error_info<struct invalid_entry_info, std::string> invalid_entry_error;
#endif // UPARSER_H_INCLUDED
