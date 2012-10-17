#ifndef UPARSERSAM_H_INCLUDED
#define UPARSERSAM_H_INCLUDED

#include "uParserBase.h"
#include <iostream>
namespace NGS
{

class uParserSam: public uParserBase
{
    class samInformation
    {
    public:
        ~samInformation() {};

    private :


    };

public :
    uParserSam();
    ~uParserSam();
    void init(const std::string& filename, bool header = false);
    void init(std::iostream* stream, bool header = false);
    void init(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
    void init(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');

    uToken getNextEntry();

private:
    static DerivedParserRegister<uParserSam> reg;
    samInformation m_Info;
    void _processSamline(std::stringstream & curSStream);
    void _parseHeader();
};

}

#endif // UPARSERSAM_H_INCLUDED
