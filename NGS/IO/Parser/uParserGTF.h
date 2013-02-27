#ifndef UPARSERGTF_H_INCLUDED
#define UPARSERGTF_H_INCLUDED

#include "uParserBase.h"
#include <iostream>
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
    static uParserBase * Create() { return new uParserGTF(); }
private:
    void _parseHeader();
    uToken _getTokenInfoFromGTFString(const std::string & line);
    boost::xpressive::sregex GTFRegex;
    const std::string GTFregString="^(\\.|[\\w_-]+)\t(\\.|[\\w_-]+)\t(\\.|[\\w_-]+)\t(\\d+)\t(\\d+)\t([-+]?[0-9]*\\.?[0-9]+|.)\t(\\+|\\-|\\.)\t([012\\.])\tgene_id\\s\"([^\"]*)\";\\stranscript_id\\s\"([^\"]*)\";.*";
};

} // End of namespace NGS


#endif // UPARSERGTF_H_INCLUDED
