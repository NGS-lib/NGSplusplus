#include "uGeneException.h"
#include "iParse.h"


    parserBase::parserBase(const std::string& filename, const std::string & type, bool header = false)
    {




    };
	parserBase::parserBase(std::iostream* stream, const std::string & type, bool header = false){



	};
	parserBase::parserBase(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t'){




	};
	parserBase::parserBase(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t'){



	};
	parserBase::~parserBase(){};

	uToken parserBase::getNextEntry(
                                return  *m_pParserBase.getNextEntry();
                                 );

