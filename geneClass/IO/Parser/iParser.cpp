#include "../../uGeneException.h"
#include "iParser.h"


    Parser::Parser(const std::string& filename, const std::string & type, bool header = false)
    {

        parserBaseFactory myFact;
        m_pParserBase=myFact.createInstance(type);
        m_pParserBase.init(filename, header);

    };
	Parser::Parser(std::iostream* stream, const std::string & type, bool header = false){

        parserBaseFactory myFact;
        m_pParserBase=myFact.createInstance(type);
        m_pParserBase.init(stream, header);

	};
	Parser::Parser(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t'){

        parserBaseFactory myFact;
        m_pParserBase=myFact.createInstance("CUSTOM");
        m_pParserBase.init(filename, fieldsNames,header,delimiter);



	};
	Parser::Parser(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t'){

        parserBaseFactory myFact;
        m_pParserBase=myFact.createInstance("CUSTOM");
        m_pParserBase.init(stream, fieldsNames,header,delimiter);

	};
	Parser::~Parser(){};

	uToken Parser::getNextEntry(
                                return  *m_pParserBase.getNextEntry();
                                 );

