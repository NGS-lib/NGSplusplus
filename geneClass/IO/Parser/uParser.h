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
#include "../uToken.h"
#include "../../uGeneException.h"
#include  "uParserBase.h"
#include "../uHeader.h"
namespace NGS {
class uHeader;
class uParser{

public :

    uParser(const std::string& filename, const std::string & type, bool header = false);
    uParser(std::iostream* stream, const std::string & type, bool header = false);
    uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    uParser(std::iostream* stream, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    ~uParser();

    bool eof() const ;
    uToken getNextEntry();

    /** \brief Get an unformated version of header (i.e.: a single string containing the whole header)
    */
    std::string getUnformatedHeader() const { return m_pParserBase->getUnformatedHeader(); }
    /** \brief Get a specific data from header.
      */
    std::string getHeaderParam(header_param name) const { return m_pParserBase->getHeaderParam(name); }
    std::vector<std::string> getHeaderParamVector(header_param name) const { return m_pParserBase->getHeaderParamVector(name); }
    /** \brief Check if there is a value associated with a given param.
      * \param header_param& name: name of the param to check.
      */
    bool isHeaderParamSet(const header_param& name) const { return m_pParserBase->isHeaderParamSet(name); }
private:
    std::shared_ptr<uParserBase> m_pParserBase=nullptr;

};
}
#endif
