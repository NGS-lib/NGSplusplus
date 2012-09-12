#ifndef IPARSER_H_INCLUDED
#define IPARSER_H_INCLUDED

#include <map>
#include <string>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "uToken.h"
#include "../uGeneException.h"
#include "uHeader.h"

class Parser{

public :

    parserBase(const std::string& filename, const std::string & type, bool header = false);
	parserBase(std::iostream* stream, const std::string & type, bool header = false);
	parserBase(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	parserBase(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	~parserBase();

    bool eof() const { return m_pIostream->peek() == EOF; }
	uToken getNextEntry();

	/** \brief Get a specific data from header.
	  */
	std::string getHeaderParam(header_param name) const { return m_headerData.getParam(name); }
	/** \brief Check if there is a value associated with a given param.
	  * \param header_param& name: name of the param to check.
	  */
	bool isHeaderParamSet(const header_param& name) const { return m_headerData.isParamSet(name); }
private:

    uHeader m_headerData
    unique_ptr<parserBase> m_pParserBase;

}


