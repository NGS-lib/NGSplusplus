#ifndef UPARSERCUSTOM_H_INCLUDED
#define UPARSERCUSTOM_H_INCLUDED

#include "uParserBed.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
#include <iostream>
//#include "uParserFactory.h"
namespace NGS
{

class uParserCustom : public uParserBed
{
public :
    uParserCustom();
    ~uParserCustom();

    virtual void init(const std::string& filename, bool header = false);
    virtual void init(std::istream* stream, bool header = false);
    void init(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    void init(std::istream* stream, const std::vector<std::string>& fieldsNames, char delimiter = '\t');
    static uParserBase * Create() { return new uParserCustom(); }
    virtual uToken getNextEntry();

private:
    std::vector<std::string> m_customFieldNames{};
    void _customParserValidateFields(const std::vector<std::string>& fieldsNames) const;
    void _customParserCopyFields(const std::vector<std::string>& fieldsNames);
    void _convertLineToTokenInfosCustom(char* line, std::stringstream& token_infos);
    bool _paramExists(const std::string& name, const std::vector<std::string>& list) const;

};

} // End of namespace NGS
#endif // UPARSERCUSTOM_H_INCLUDED
