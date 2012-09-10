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
	uParser(const std::string& filename, file_type type, bool header = false);
	uParser(std::iostream* stream, file_type type, bool header = false);
	uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	uParser(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	~uParser();
	uParser& operator=(const uParser& copyFrom) = delete;
	uParser(const uParser&) = delete;
	bool eof() const { return m_pIostream->peek() == EOF; }
	uToken getNextEntry();

private:
	std::iostream* m_pIostream = nullptr;
	file_type m_fileType = file_type::AUTO;
	char m_delimiter = '\t';
	bool m_header = false;
//	bool m_firstToken = true;
	std::vector<std::string> m_customFieldNames{};
	//TODO: To avoid using delete on m_pIostream if ifstream constructor was used. Is there a better way?
	bool m_dynamicStream = false;

	void _fetchHeader();
	void _fetchUnspecifiedHeader();
	void _pushBackLine(char* line);
	uToken _getNextEntryBed();
	void _convertLineToTokenInfosBed(char* line, std::stringstream& token_infos);
	uToken _getNextEntrySam();
	uToken _getNextEntryCustom();
	void _convertLineToTokenInfosCustom(char* line, std::stringstream& token_infos);
	bool _paramExists(const std::string& name, const std::vector<std::string>& list) const;
	void _customParserValidateFields(const std::vector<std::string>& fieldNames);
	void _customParserCopyFields(const std::vector<std::string>& fieldsNames);
};

#endif // UPARSER_H_INCLUDED
