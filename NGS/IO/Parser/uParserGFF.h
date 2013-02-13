#ifndef UPARSERGFF_H_INCLUDED
#define UPARSERGFF_H_INCLUDED

#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>
#include "../../boost-include/boost/xpressive/xpressive.hpp"
#include "uParserBase.h"
namespace NGS
{

    class uParserGFF : public uParserBase
    {
    public :
        uParserGFF();
        ~uParserGFF();

        virtual void init(const std::string& filename, bool header = false);
        virtual void init(std::istream* stream, bool header = false);
        uToken getNextEntry();
        static uParserBase * Create() { return new uParserGFF(); }
    private:
        uToken _getTokenFromGFFString(const std::string & line);
      //  static DerivedParserRegister<uParserGFF> reg;
        boost::xpressive::sregex GFFRegex;
        const std::string GFFregString="^(\\.|[\\w_-]+)\t(\\.|[\\w_-]+)\t(\\.|[\\w_-]+)\t(\\d+)\t(\\d+)\t([-+]?[0-9]*\\.?[0-9]+|.)\t(\\+|\\-|\\.)\t([012\\.])(?:\t(.+))?";
    };

} // End of namespace NGS

#endif // UPARSERGFF_H_INCLUDED
