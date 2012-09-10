#ifndef UPARSER_H_INCLUDED
#define UPARSER_H_INCLUDED

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "uToken.h"
#include "uGeneException.h"

enum class file_type{ BED, SAM, CUSTOM ,WIG};

class typeInformation{

public:
  virtual ~typeInformation() {};
};

class samInformation : public typeInformation{
public:
    ~samInformation() {};
};
class wigInformation : public typeInformation{
public:

     ~wigInformation() {};
    enum class stepType{ FIXED,VARIABLE, NA };



    stepType getStep() {return m_step;};
    void setStep(stepType p_step){m_step=p_step;};

    std::string getChrom(){return m_chrom;};
    void setChrom(std::string p_chrom){m_chrom=p_chrom;};
    void setSpan(long int p_span){m_span=p_span;};
    long int getSpan(){return m_span ;};
    private :

    stepType m_step=stepType::NA;
    std::string m_chrom="";
    long int m_span=0;
};
class bedInformation : public typeInformation{
};
class customInformation : public typeInformation{
};

//std::unique_ptr<typeInformation

class uParser {
public:
	uParser(const std::string& filename, file_type type, bool header = false);
	uParser(std::iostream* stream, file_type type, bool header = false);
	uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	uParser(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	~uParser();
	uParser& operator=(const uParser& copyFrom) = delete;
	uParser(const uParser&) = delete;
	bool eof() const { return m_pIostream->peek() == EOF; }
	uToken getNextEntry();

private:
	std::istream* m_pIostream = nullptr;
	file_type m_fileType;

	char m_delimiter = '\t';
	bool m_header = false;
//	bool m_firstToken = true;
	std::vector<std::string> m_customFieldNames{};
	//TODO: To avoid using delete on m_pIostream if ifstream constructor was used. Is there a better way?
	bool m_dynamicStream = false;

    std::unique_ptr<typeInformation> m_info=nullptr;
    std::unique_ptr<typeInformation> _makeTypeInfo(file_type type);


    void _processWigDeclaration(std::stringstream & curSStream);


	void _fetchHeader();
	void _fetchUnspecifiedHeader();
	void _pushBackLine(char* line);
	uToken _getNextEntryBed();
	void _convertLineToTokenInfosBed(char* line, std::stringstream& token_infos);
	uToken _getNextEntrySam();
	uToken _getNextEntryCustom();
    uToken _getNextEntryWig();
	void _convertLineToTokenInfosCustom(char* line, std::stringstream& token_infos);

	bool _paramExists(const std::string& name, const std::vector<std::string>& list) const;
	void _customParserValidateFields(const std::vector<std::string>& fieldNames);
	void _customParserCopyFields(const std::vector<std::string>& fieldsNames);
};

#endif // UPARSER_H_INCLUDED
