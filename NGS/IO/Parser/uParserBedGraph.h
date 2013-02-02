#ifndef UPARSERBEDGRAPH_INCLUDED
#define UPARSERBEDGRAPH_INCLUDED

#include "uParserBase.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>
#include "../../boost-include/boost/xpressive/xpressive.hpp"
//#include "uParserFactory.h"
namespace NGS
{

class uParserBedGraph : public uParserBase
{
public :
    uParserBedGraph();
    ~uParserBedGraph();

    virtual void init(const std::string& filename, bool header = true);
    virtual void init(std::istream* stream, bool header = true);
    uToken getNextEntry();
     static uParserBase * Create() { return new uParserBedGraph(); }
private:
    uToken _getTokenFromBedGraphString(const std::string & line);
    bool _parseHeader();
//    static DerivedParserRegister<uParserBedGraph> reg;
    bool m_headerFound=false;
	std::stringstream m_hBuffer;
	const std::string s_bedGraphHeader="track type=bedGraph";
    std::vector<std::string> m_tokens;
};

} // End of namespace NGS

#endif // UPARSERBEDGRAPH_INCLUDED
