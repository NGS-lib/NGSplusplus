#include "uParserBase.h"

namespace NGS 
{
std::map<std::string, std::function<uParserBase*()> > *uParserBaseFactory::mapItem;

uParserBase::uParserBase(){};

bool uParserBase::eof() const
{ 
    return m_pIostream->peek() == EOF; 
}

}
