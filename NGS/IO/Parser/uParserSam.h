#ifndef UPARSERSAM_H_INCLUDED
#define UPARSERSAM_H_INCLUDED

#include "uParserBase.h"
#include "../../utility/utility.h"
#include <iostream>
#include "../../boost-include/boost/xpressive/xpressive.hpp"
//#include "uParserFactory.h"
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
	static uParserBase * Create() { return new uParserSam(); }
private:
 //   static DerivedParserRegister<uParserSam> reg;
    samInformation m_Info;
    void _parseHeader();
    /**< String for dynamic parsing of Sam */
	const std::string SamRegString="^([!-?A-~]{1,255})\t([\\d]+)\t(\\*|[!-()+-<>-~][!-~]*)\t(\\d+)\t(\\d+)\t(\\*|(?:[0-9]+[MIDNSHPX=])+)\t(\\*|=|[!-()+-<>-~][!-~]*)\t(\\d+)\t(-{0,1}\\d+)\t(\\*|[A-Za-z=.]+)\t([!-~]+)(?:\t(.+))?";

	/**< Necessary for static regex */
	boost::xpressive::smatch what;
	boost::xpressive::mark_tag s10;
	boost::xpressive::mark_tag s11;

	boost::xpressive::mark_tag s12;
	boost::xpressive::sregex staticSam= boost::xpressive::bos >>(boost::xpressive::s1 =(boost::xpressive::repeat<1,255>(boost::xpressive::set[ boost::xpressive::range('!','?') | boost::xpressive::range('A','~') ]) ) )>>"\t">>(boost::xpressive::s2 =(+boost::xpressive::_d ))>>"\t">>
	(boost::xpressive::s3=('*'|boost::xpressive::set[boost::xpressive::range('!','(')|')'|boost::xpressive::range('+','<')|boost::xpressive::range('>','~')]>>*boost::xpressive::set[boost::xpressive::range('!','~')] ))
	>>"\t">>(boost::xpressive::s4 =(+boost::xpressive::_d ))>>"\t">>(boost::xpressive::s5 =(+boost::xpressive::_d ))>>"\t">>
	(boost::xpressive::s6=('*'|+(+boost::xpressive::set[boost::xpressive::range('0','9')]>>+(boost::xpressive::set='M','I','D','N','S','H','P','X','='))))
	>>"\t">>(boost::xpressive::s7=((boost::xpressive::set= '*','=') | (boost::xpressive::set[boost::xpressive::range('!','(')|')'|boost::xpressive::range('+','<')|boost::xpressive::range('>','~')]>>*boost::xpressive::set[boost::xpressive::range('!','~')])  ))
	>>"\t">>(boost::xpressive::s8 =(+boost::xpressive::_d ))>>"\t">>(boost::xpressive::s9 =(boost::xpressive::repeat<0,1>('-')>>+boost::xpressive::_d))>>"\t">>(s10=('*'|+boost::xpressive::set[boost::xpressive::range('A','Z')|boost::xpressive::range('a','z')|(boost::xpressive::set='=','.')]) )
	>>"\t">>(s11=+boost::xpressive::set[boost::xpressive::range('!','~')])>>!(("\t")>>(s12=(+boost::xpressive::_)));
};

}

#endif // UPARSERSAM_H_INCLUDED
