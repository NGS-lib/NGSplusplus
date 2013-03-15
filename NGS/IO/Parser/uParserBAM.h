#ifndef UPARSEBAM_H_INCLUDED
#define UPARSEBAM_H_INCLUDED
#include "uParserBase.h"
#include <iostream>
#include "api/BamReader.h"
namespace NGS
{

class uParserBAM: public uParserBase
{

public :
    uParserBAM();
    ~uParserBAM();
    void init(const std::string& filename, bool header = false);
    void init(std::istream* stream, bool header = false);

    uToken getNextEntry();
	static uParserBase * Create() { return new uParserBAM(); }
	bool eof();
private:
 //   static DerivedParserRegister<uParserBAM> reg;
    BamTools::BamReader m_BamReader;
    bool m_IsBuffer=false;
    BamTools::BamAlignment m_BufferAlignement;
    void _parseHeader();
    /**< String for dynamic parsing of Sam */
};

}


#endif // UPARSEBAM_H_INCLUDED
