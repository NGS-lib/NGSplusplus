#ifndef UPARSERSAM_H_INCLUDED
#define UPARSERSAM_H_INCLUDED

#include "uParserBase.h"
#include <iostream>
#include "../../boost-include/boost/xpressive/xpressive.hpp"
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
    void init(std::istream* stream, bool header = false);

    uToken getNextEntry();
	uToken getNextEntryWithRegex();
private:
    static DerivedParserRegister<uParserSam> reg;
    samInformation m_Info;
	boost::xpressive::sregex SAMRegex;
    void _parseHeader();
	const std::string SamRegString="^([!-?A-~]{1,255})\t([\\d]+)\t(\\*|[!-()+-<>-~][!-~]*)\t(\\d+)\t(\\d+)\t(\\*|(?:[0-9]+[MIDNSHPX=])+)\t(\\*|=|[!-()+-<>-~][!-~]*)\t(\\d+)\t(-{0,1}\\d+)\t(\\*|[A-Za-z=.]+)\t([!-~]+)(?:\t(.+))?";

};

}

#endif // UPARSERSAM_H_INCLUDED
