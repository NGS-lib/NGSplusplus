#include "../../uGeneException.h"
#include "uParser.h"
#include "uParserBase.h"
namespace NGS {

    uParser::uParser(const std::string& filename, const std::string & type, bool header)
    {

        uParserBaseFactory myFact;
        m_pParserBase=myFact.createInstance(type);
        m_pParserBase->init(filename, header);

    };
	uParser::uParser(std::iostream* stream, const std::string & type, bool header){

        uParserBaseFactory myFact;
        m_pParserBase=myFact.createInstance(type);
        m_pParserBase->init(stream, header);

	};
	uParser::uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header, char delimiter){

        uParserBaseFactory myFact;
        m_pParserBase=myFact.createInstance("CUSTOM");
        m_pParserBase->init(filename, fieldsNames,header,delimiter);

	};
	uParser::uParser(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header, char delimiter){

        uParserBaseFactory myFact;
        m_pParserBase=myFact.createInstance("CUSTOM");
        m_pParserBase->init(stream, fieldsNames,header,delimiter);

	};
	uParser::~uParser(){};

     bool uParser::eof() const {
         return m_pParserBase->eof();
         };

    uToken uParser::getNextEntry(){
	 return  m_pParserBase->getNextEntry();
	};
}
