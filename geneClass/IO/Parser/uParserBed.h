#ifndef UPARSERBED_H_INCLUDED
#define UPARSERBED_H_INCLUDED

#include "uParserBase.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>

namespace NGS
{

class uParserBed : public uParserBase
{
public :
    uParserBed();
    ~uParserBed();

    void init(const std::string& filename, bool header = false);
    void init(std::iostream* stream, bool header = false);
    void init(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    void init(std::iostream* stream, const std::vector<std::string>& fieldsNames, char delimiter = '\t');

    uToken getNextEntry();

protected:
    char m_delimiter = '\t';

private:
    void _parseHeader();
    int _countColumns(char* line) const;
    void _validateColumnNumber(int numberOfColumn) const;
    void _convertLineToTokenInfosBed(char* line, std::stringstream& token_infos);
    void _pushBackLine(char* line);
    std::string _getNextEntry(char* line);

    int m_numberOfColumn = 0;
    bool m_headerParsed = false;
    static DerivedParserRegister<uParserBed> reg;
};

} // End of namespace NGS
#endif // UPARSERBED_H_INCLUDED
