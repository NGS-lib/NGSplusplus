#ifndef UPARSER_H_INCLUDED
#define UPARSER_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include "uToken.h"

enum class file_type{ AUTO, BED, SAM };

class uParser {
public:
	uParser(const std::string& filename, file_type type);
	uParser(std::istream* stream, file_type type);
	~uParser();
	uParser& operator=(const uParser& copFrom)=delete;
	uParser(const uParser&) = delete;

	bool eof() const { return m_pIstream->peek() == EOF; }
	uToken getNextEntry();

private:
//	std::istream m_istream;
//	std::istream* m_ifstream;
	std::istream* m_pIstream=nullptr;
	file_type m_fileType=file_type::AUTO;
	uToken _getNextEntryBed();
	uToken _getNextEntrySam();
	//TODO: To avoid using delete on m_pIstream if ifstream constructor was used. Is there a better way?
	bool m_dynamicStream=false;
};

/**<  uParser exceptions */
struct uParser_exception_base : virtual std::exception, virtual boost::exception {};

struct end_of_file_throw : virtual uToken_exception_base{};
typedef boost::error_info<struct end_of_file_info, std::string> end_of_file_error;
#endif // UPARSER_H_INCLUDED
