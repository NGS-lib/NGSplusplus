#ifndef UPARSERWIG_H_INCLUDED
#define UPARSERWIG_H_INCLUDED

#include "iParserBase.h"
#include <iostream>

namespace NGS{

class uParserWig: public uParserBase
{
    public :
    uParserWig();
    ~uParserWig(){};
    void init(const std::string& filename, bool header = false);
	void init(std::iostream* stream, bool header = false);
	void init(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	void init(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');

	uToken getNextEntry();

private:
    static DerivedRegister<uParserWig> reg;
};

}
#endif // UPARSERWIG_H_INCLUDED
