#include "uParserWig.h"

namespace NGS{
uParserWig::uParserWig () : uParserBase()
{


}

    void uParserWig::init(const std::string& filename, bool header )
    {



    }
	void uParserWig::init(std::iostream* stream, bool header )
	{


	}
	void uParserWig::init(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header , char delimiter)
	{


	}
	void uParserWig::init(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header , char delimiter )
	{


	}

    uToken uParserWig::getNextEntry(){


    }

DerivedRegister<uParserWig> uParserWig::reg("WIG");

}
