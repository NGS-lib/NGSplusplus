#ifndef UPARSER_H_INCLUDED
#define UPARSER_H_INCLUDED

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "uToken.h"
#include "uGeneException.h"

enum class file_type{ AUTO, BED, SAM, CUSTOM };

class uParser {
public:
	uParser(const std::string& filename, file_type type);
	uParser(std::istream* stream, file_type type);
	uParser(const std::vector<std::string>& columnNames, char delimiter = '\t');
	~uParser();
	uParser& operator=(const uParser& copyFrom) = delete;
	uParser(const uParser&) = delete;
	bool eof() const { return m_pIstream->peek() == EOF; }
	uToken getNextEntry();

private:
	std::istream* m_pIstream = nullptr;
	file_type m_fileType = file_type::AUTO;
	char m_delimiter = '\t';
	std::vector<std::string> m_customFieldNames;
	//TODO: To avoid using delete on m_pIstream if ifstream constructor was used. Is there a better way?
	bool m_dynamicStream = false;
	
	uToken _getNextEntryBed();
	uToken _getNextEntrySam();
	void _customParserValidateFields(const std::vector<std::string>& fieldNames);
	void _customParserCopyFields(const std::vector<std::string>& fieldsNames);
};

#endif // UPARSER_H_INCLUDED
