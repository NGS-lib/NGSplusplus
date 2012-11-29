#ifndef UPARSERGTF_H_INCLUDED
#define UPARSERGTF_H_INCLUDED

#include "uParserBase.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>
#include "../../boost-include/boost/xpressive/xpressive.hpp"
namespace NGS
{

class uParserGTF : public uParserBase
{
public :
    uParserGTF();
    ~uParserGTF();

    virtual void init(const std::string& filename, bool header = false);
    virtual void init(std::istream* stream, bool header = false);
    uToken getNextEntry();

private:
    void _getTokenInfoFromGTFString(std::string line, std::stringstream& token_infos);
    static DerivedParserRegister<uParserGTF> reg;
    boost::xpressive::sregex GTFRegex;
    const std::string GTFregString="(\\w+)\t(\\w+)\t(\\w+)\t(\\d+)\t(\\d+)\t([-+]?[0-9]*\\.?[0-9]+|.)\t(\\+|\\-|\\.)\t([012])(?:\t(.+))?";
};

} // End of namespace NGS


#endif // UPARSERGTF_H_INCLUDED
