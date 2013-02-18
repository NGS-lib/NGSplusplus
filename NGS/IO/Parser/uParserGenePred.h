#ifndef _H_INCLUDED
#define _H_INCLUDED
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>
//#include "../../boost-include/boost/xpressive/xpressive.hpp"
#include "uParserBase.h"
namespace NGS
{

    class uParserGenePred : public uParserBase
    {
    public :
        uParserGenePred();
        ~uParserGenePred();

        virtual void init(const std::string& filename, bool header = false);
        virtual void init(std::istream* stream, bool header = false);
        uToken getNextEntry();
        static uParserBase * Create() { return new uParserGenePred(); }
    private:
        uToken _getTokenFromGenePredString(const std::string & line);
        void  _parseHeader();
        boost::xpressive::sregex UCSCGenePredRegex;
        boost::xpressive::sregex UCSCGenePredRegexInit;
       // const std::string GenePredregStringPriortest="(.*)";
        const std::string GenePredregStringPrior="^([.\\w_-]+)\t([.\\w_-]+)\t(\\+|\\-)\t(\\d+)\t(\\d+)\t(\\d+)\t(\\d+)\t(\\d+)\t(.*)";
        const std::string GenePredregStringPart1="^([.\\w_-]+)\t([.\\w_-]+)\t(\\+|\\-)\t(\\d+)\t(\\d+)\t(\\d+)\t(\\d+)\t(\\d+)";
        const std::string GenePredRegExon ="(\\d+),";

    };

} // End of namespace NGS

#endif // _H_INCLUDED
