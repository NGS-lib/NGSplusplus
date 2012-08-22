#ifndef UPARSER_H_INCLUDED
#define UPARSER_H_INCLUDED

#include <string>
#include <fstream>
#include "uToken.h"

enum class file_type{ BED, SAM };

class uParser {
public:
	uParser(const std::string& filename, file_type type);
	uToken getNextEntry();

private:
	std::ifstream m_ifstream;
	file_type m_fileType;
	uToken _getNextEntryBed();
};

/**<  uParser exceptions */
struct uParser_exception_base : virtual std::exception, virtual boost::exception {};

struct end_of_file_throw : virtual uToken_exception_base{};
typedef boost::error_info<struct end_of_file_info, std::string> end_of_file_error;
#endif // UPARSER_H_INCLUDED
