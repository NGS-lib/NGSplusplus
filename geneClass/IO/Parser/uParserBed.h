#ifndef UPARSERBED_H_INCLUDED
#define UPARSERBED_H_INCLUDED

#include "iParserBase.h"
#include <iostream>

namespace NGS
{

class uParserBed: public uParserBase
{
    class wigInformation
    {
    public:

    private :
        long int m_curPos=0;
        stepType m_stepType=stepType::NA;
        std::string m_chrom="";
        long int m_span=-1;
    };

public :
    uParserBed();
    ~uParserBed();
    void init(const std::string& filename, bool header = false);
    void init(std::iostream* stream, bool header = false);
    void init(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
    void init(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');

    uToken getNextEntry();

private:
    static DerivedParserRegister<uParserBed> reg;
   // wigInformation m_Info;
    void _processFixedWigLine(std::stringstream & curSStream);
    void _processVariabledWigLine(std::stringstream & curSStream);
};

}


#endif // UPARSERBED_H_INCLUDED
