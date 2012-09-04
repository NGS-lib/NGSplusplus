#ifndef UPARSER_H_INCLUDED
#define UPARSER_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include "uToken.h"
#include "uGeneException.h"

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
	std::istream* m_pIstream=nullptr;
	file_type m_fileType=file_type::AUTO;
	uToken _getNextEntryBed();
	uToken _getNextEntrySam();
	//TODO: To avoid using delete on m_pIstream if ifstream constructor was used. Is there a better way?
	bool m_dynamicStream=false;
};

#endif // UPARSER_H_INCLUDED
