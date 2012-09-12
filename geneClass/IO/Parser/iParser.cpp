#include "../../uGeneException.h"
#include "iParser.h"
#include "iParserBase.h"
namespace NGS {

    Parser::Parser(const std::string& filename, const std::string & type, bool header)
    {

        uParserBaseFactory myFact;
        m_pParserBase=myFact.createInstance(type);
        m_pParserBase->init(filename, header);

    };
	Parser::Parser(std::iostream* stream, const std::string & type, bool header){

        uParserBaseFactory myFact;
        m_pParserBase=myFact.createInstance(type);
        m_pParserBase->init(stream, header);

	};
	Parser::Parser(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header, char delimiter){

        uParserBaseFactory myFact;
        m_pParserBase=myFact.createInstance("CUSTOM");
        m_pParserBase->init(filename, fieldsNames,header,delimiter);

	};
	Parser::Parser(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header, char delimiter){

        uParserBaseFactory myFact;
        m_pParserBase=myFact.createInstance("CUSTOM");
        m_pParserBase->init(stream, fieldsNames,header,delimiter);

	};
	Parser::~Parser(){};

     bool Parser::eof() const {
         return m_pParserBase->eof();
         };

    uToken Parser::getNextEntry(){
	 return  m_pParserBase->getNextEntry();
	};
}
