#include "../../uGeneException.h"
#include "uParser.h"
#include "uParserBase.h"
namespace NGS
{

uParser::uParser(const std::string& filename, const std::string & type, bool header)
{

    uParserBaseFactory myFact;
    try
    {
        m_pParserBase=myFact.createInstance(type);
        m_pParserBase->init(filename, header);

    }
    catch(...)
    {
        throw;
        //  std::cerr <<fetchStringError(e);
    }
};
uParser::uParser(std::iostream* stream, const std::string & type, bool header)
{

    uParserBaseFactory myFact;
    m_pParserBase=myFact.createInstance(type);
    m_pParserBase->init(stream, header);

};
uParser::uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter)
{

    uParserBaseFactory myFact;
    m_pParserBase=myFact.createInstance("CUSTOM");
    m_pParserBase->init(filename, fieldsNames, delimiter);

};
uParser::uParser(std::iostream* stream, const std::vector<std::string>& fieldsNames, char delimiter)
{

    uParserBaseFactory myFact;
    m_pParserBase=myFact.createInstance("CUSTOM");
    m_pParserBase->init(stream, fieldsNames, delimiter);

};
uParser::~uParser() {};

bool uParser::eof() const
{
    return m_pParserBase->eof();
};

uToken uParser::getNextEntry()
{
    return  m_pParserBase->getNextEntry();
};
}
