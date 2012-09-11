#include "uParser.h"
#include "uGeneException.h"
/** \brief uParser constructor with filename
 * \param std::string filename: Name of the file to load
 * \param file_type type: Currently supported formats: BED, SAM
 */
uParser::uParser(const std::string& filename, file_type type, bool header):m_fileType(type),m_header(header)
{
    std::ifstream* ifs = new std::ifstream(filename.c_str(), std::ifstream::in);
    if (!ifs->is_open())
    {
        std::string error = "Error opening file: " + filename;
        throw std::runtime_error(error.c_str());
    }
    else
    {
        m_pIostream = ifs;
    }

    m_info=_makeTypeInfo(type);
    m_header = header;
    if (m_header == true)
    {
        _fetchHeader();
    }
    m_dynamicStream = true;
}


/** \brief uParser constructor with istream
 * \param std::istream* stream: stream to load the data from
 * \param file_type type: Currently supported formats: BED, SAM
 */
uParser::uParser(std::iostream* stream, file_type type, bool header):m_fileType(type),m_header(header)
{
    m_pIostream = stream;
    m_info=_makeTypeInfo(type);
    m_header = header;
    if (m_header == true)
    {
        _fetchHeader();
    }
    m_dynamicStream = false;
}

/** \brief uParser constructor with custom file type
 * \param const std::vector<string> columnNames: The name of every field in the custom format, in the SAME ORDER as they appear in the file. Must be of the string version of token_param (see uToken.h). Mandatory fields are: CHR, START_POS and END_POS.
 * \param char delimiter: The delimiter between each field.
 */
uParser::uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header, char delimiter):m_fileType(file_type::CUSTOM)
{
    /**< Check if filename is valid, then open it */
    std::fstream* ifs = new std::fstream(filename.c_str(), std::fstream::in);
    if (!ifs->is_open())
    {
        std::string error = "Error opening file: " + filename;
        throw std::runtime_error(error.c_str());
    }
    else
    {
        m_pIostream = ifs;
    }
    /**< Check if fields are in a valid format */
    try
    {
        _customParserValidateFields(fieldsNames);
        _customParserCopyFields(fieldsNames);
    }
    catch(customParser_missing_mandatory_values& e)
    {
        throw e;
    }
    /**< Set other parameters */

    m_delimiter = delimiter;
    m_header = header;
    if (m_header == true)
    {
        _fetchHeader();
    }
    m_dynamicStream = true;
}

uParser::uParser(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header, char delimiter):m_fileType(file_type::CUSTOM)
{

    /**< Check if fields are in a valid format */
    try
    {
        _customParserValidateFields(fieldsNames);
        _customParserCopyFields(fieldsNames);
    }
    catch(customParser_missing_mandatory_values& e)
    {
        throw e;
    }
    /**< Set other parameters */
    m_pIostream = stream;
    m_delimiter = delimiter;
    m_header = header;
    if (m_header == true)
    {
        _fetchHeader();
    }
    m_dynamicStream = false;
}

/** \brief uParser destructor
 */
uParser::~uParser()
{
    if (m_dynamicStream == true)
    {
        delete m_pIostream;
    }
    m_pIostream = NULL;
}

/** \brief Generic function to fetch an entry from file, will adapt itself to file_type given in constructor
 */
uToken uParser::getNextEntry()
{
    try
    {
        switch(m_fileType)
        {
        case file_type::BED:
            try
            {
                uToken token = _getNextEntryBed();
                return token;
            }
            catch(invalid_uToken_throw& e)
            {
#ifdef DEBUG
                std::string trace = fetchStringError(e);
                std::cerr << "Invalid uToken: " << trace << std::endl;
#endif
                throw e;
            }
            catch(end_of_file_throw& e)
            {
                throw e;
            }
            break;

        case file_type::SAM:
            try
            {
                uToken token = _getNextEntrySam();
                return token;
            }
            catch(invalid_uToken_throw& e)
            {
#ifdef DEBUG
                std::string trace = fetchStringError(e);
                std::cerr << "Invalid uToken: " << trace << std::endl;
#endif
                throw e;
            }
            break;

        case file_type::CUSTOM:
            try
            {
                uToken token = _getNextEntryCustom();
                return token;
            }
            catch(invalid_uToken_throw& e)
            {
#ifdef DEBUG
                std::string trace = fetchStringError(e);
                std::cerr << "Invalid uToken: " << trace << std::endl;
#endif
                throw e;
            }
        case file_type::WIG:
            try
            {
                uToken token = _getNextEntryWig();
                return token;
            }
            catch(invalid_uToken_throw& e)
            {
#ifdef DEBUG
                std::string trace = fetchStringError(e);
                std::cerr << "Invalid uToken: " << trace << std::endl;
#endif
                throw e;
            }
        default:
            throw uParser_exception_base()<< string_error("Invalid fileType in getNextEntry case");
            break;
        }
    }
    catch(invalid_uToken_throw &e)
    {
        throw e;
    }
}

// TODO: class uHeader to return from this function. Derived from uToken?
/** \brief Fetch header differently based on file type
 */
void uParser::_fetchHeader()
{
    switch(m_fileType)
    {
    case file_type::BED:
    case file_type::CUSTOM:
        _fetchUnspecifiedHeader();
        break;
        // TODO: Write header parser for SAM format
    case file_type::SAM:
        break;
    case file_type::WIG:
        _fetchWigHeader();
        break;

    default:
        break;
    }
}

/** \brief return the appropriate typeInformation to store when parsing a filetype
 * \return unique_ptr object, holding pointer to appropriate object
 */
std::shared_ptr<typeInformation>  uParser::_makeTypeInfo(file_type type)
{

    switch(m_fileType)
    {
    case file_type::BED:
        return std::shared_ptr<typeInformation>(new bedInformation);
    case file_type::CUSTOM:
        return std::shared_ptr<typeInformation>(new customInformation);
    case file_type::SAM:
        return std::shared_ptr<typeInformation>(new samInformation);
    case file_type::WIG:
        return std::shared_ptr<typeInformation>(new wigInformation);
    }


    throw uParser_exception_base()<<string_error("Throwing in _makeTypeInfo, reached non-valid area in function \n");

}

/** \brief Simply fetch header without parsing it
 */
void uParser::_fetchUnspecifiedHeader()
{
    bool headerFetched = false;
    std::string unformatedHeader;
    /**< We parse a line at a time until we get a valid token, then we return the token back in the stream */
    do
    {
        char line[4096];
        if (m_pIostream->getline(line, 4096))
        {
            std::stringstream token_infos;
            if (m_fileType == file_type::BED)
            {
                _convertLineToTokenInfosBed(line, token_infos);
            }
            else if (m_fileType == file_type::CUSTOM)
            {
                _convertLineToTokenInfosCustom(line, token_infos);
            }
            try
            {
                uToken token(token_infos);
                /**< When we have a valid token, we consider that all the header have beed fetched */
                /**< We put token back into stream, and exit function */
                _pushBackLine(line);
                headerFetched = true;
            }
            /**< If there is a throw, we are not at the first valid token line, so we add line to header object */
            catch(invalid_uToken_throw& e) { }
            catch(invalid_value_throw& e) { }
            if (headerFetched == false)
            {
                std::string s(line);
                unformatedHeader += s;
            }
        }
        else
        {
            std::cout << "Reached end of file" << std::endl;
#ifdef DEBUG
            std::cerr << "Reached end of file." << std::endl;
#endif
            end_of_file_throw e;
            e << string_error("Reached end of file.");
            throw e;
        }
    }
    while(headerFetched == false);
    m_headerData.setUnformatedHeader(unformatedHeader);
}

/** \brief Same as unspecified, but we stop at a track line. If none met, invalid file
 */
void uParser::_fetchWigHeader()
{
    bool headerFetched = false;
    std::string unformatedHeader;
    /**< We parse a line at a time until we get a valid track definition line, then we return the track line back in the stream */
    do
    {
        char line[4096];
        if (m_pIostream->getline(line, 4096))
        {
            std::stringstream ss;
            ss << line;
            std::string first_token;
            ss >> first_token;
            if((first_token=="variableStep")||(first_token=="fixedStep"))
            {
                headerFetched = true;
                _pushBackLine(line);
            }
            else
            {
                std::string s(line);
                unformatedHeader += s;
            }
        }
        else
        {
            std::cout << "Reached end of file" << std::endl;
#ifdef DEBUG
            std::cerr << "Reached end of file." << std::endl;
#endif
            end_of_file_throw e;
            e << string_error("Reached end of file.");
            throw e;
        }
    }
    while(headerFetched == false);
    m_headerData.setUnformatedHeader(unformatedHeader);
}


void uParser::_pushBackLine(char* line)
{
    m_pIostream->putback('\n');
    std::string str_line(line);
    for (int i = (int)(str_line.size()) - 1; i >= 0; i--)
    {
        m_pIostream->putback(line[i]);
    }

}

// TODO: Return pointer instead of copy of token, it could save a lot of time!
/** \brief Specific loader for BED file (See genome.ucsc.edu/FAQ/FAQformat.html#format1 for bed description)
 * \return uToken: If all the parameters in the entry are valid a uToken object is returned.
 */
uToken uParser::_getNextEntryBed()
{
    char line[4096];
    if (m_pIostream->getline(line, 4096))
    {
        /**< We start by fetching the infos from the line */
        std::stringstream token_infos;
        _convertLineToTokenInfosBed(line, token_infos);
        /**< We try to create a token with the infos that were fetched from the line */
        /**< If it doesn't work andit's the first token, we don't throw an error. Instead, we try again with another line */
        try
        {
            uToken token(token_infos);
            return token;
        }
        catch(invalid_uToken_throw& e)
        {
            throw e;
        }
        catch(invalid_value_throw& e)
        {
            throw e;
        }
    }
    else
    {
#ifdef DEBUG
        std::cerr << "Reached end of file." << std::endl;
#endif
        end_of_file_throw e;
        e << string_error("Reached end of file.");
        throw e;
    }
    std::cerr <<"Fatal error in _getNextEntryBed(), should not reach here" <<std::endl;
    abort();
}

void uParser::_convertLineToTokenInfosBed(char* line, std::stringstream& token_infos)
{
    std::stringstream ss;
    ss << line;
    std::string chr;
    std::string start_pos;
    std::string end_pos;
    std::string score;
    std::string seq_name;
    std::string strand;
    ss >> chr >> start_pos >> end_pos >> seq_name >> score >> strand;
    token_infos << "CHR\t" << chr << "\n";
    token_infos << "START_POS\t" << start_pos << "\n";
    token_infos << "END_POS\t" << end_pos << "\n";
    token_infos << "SEQ_NAME\t" << seq_name << "\n";
    /**< If there was no strand info, we don't add an empty string */
    if (strand.size() != 0)
    {
        token_infos << "STRAND\t" << strand << "\n";
    }
}

/** \brief Specific loader for SAM file (See samtools.sourceforge.net for SAM description)
 * \return uToken: If all the parameters in the entry are valid a uToken object is returned.
 */
uToken uParser::_getNextEntrySam()
{
    try
    {
        char line[4096];
        if (m_pIostream->getline(line, 4096))
        {
            std::stringstream ss;
            ss << line;

            //String to test for int
            std::string flag;
            std::string start_pos;
            std::string MAPQual;
            std::string pNext;
            std::string Tlen;

            std::string chr;

            std::string end_pos;
            std::string score;
            std::string seq_name;
            std::string qual;
            std::string seq;
            std::string cigar;
            std::string RNext;
            std::stringstream token_infos;

            ss >> seq_name >> flag >> chr >> start_pos >> MAPQual >> cigar>>RNext>>pNext>>Tlen>>seq>>qual;

            token_infos << "CHR\t" << chr << "\n";
            token_infos << "START_POS\t" << start_pos << "\n";
            //token_infos << "END_POS\t" << end_pos << "\n";
            token_infos << "FLAG\t" << flag << "\n";
            token_infos << "SEQ_NAME\t" << seq_name << "\n";
            token_infos << "MAP_SCORE\t" << MAPQual << "\n";
            token_infos << "SEQUENCE\t" << seq << "\n";
            token_infos << "CIGAR\t" << cigar << "\n";
            token_infos << "PHRED_SCORE\t" << qual << "\n";

            uToken token(token_infos);
            return token;
        }
        else
        {
#ifdef DEBUG
            std::cerr << "Reached end of file." << std::endl;
#endif
            end_of_file_throw e;
            e << string_error("Reached end of file.");
            throw e;
        }
    }
    catch(invalid_uToken_throw& e)
    {
        throw e;
    }
}

uToken uParser::_getNextEntryCustom()
{
    char line[4096];
    if (m_pIostream->getline(line, 4096))
    {
        /**< We start by fetching the infos in the line */
        std::stringstream token_infos;
        _convertLineToTokenInfosCustom(line, token_infos);
        /**< Then we try to create a token with that info */
        /**< If it doesn't work andit's the first token, we don't throw an error. Instead, we try again with another line */
        try
        {
            uToken token(token_infos);
            return token;
        }
        catch(invalid_uToken_throw& e)
        {
            throw e;
        }
        catch(invalid_value_throw& e)
        {
            throw e;
        }
    }
    else
    {
#ifdef DEBUG
        std::cerr << "Reached end of file." << std::endl;
#endif
        end_of_file_throw e;
        e << string_error("Reached end of file.");
        throw e;
    }
    std::cerr <<"Fatal error in _getNextEntryCustom(), should not reach here" <<std::endl;
    abort();
}

void uParser::_convertLineToTokenInfosCustom(char* line, std::stringstream& token_infos)
{
    char l[4096];
    for (int i = 0; i < 4096; i++)
    {
        l[i] = line[i];
    }
    char* current;
    current = strtok(l, &m_delimiter);
    for(size_t i = 0; i < m_customFieldNames.size(); i++)
    {
        if (m_customFieldNames[i] != "NA")
        {
            token_infos << m_customFieldNames[i] << "\t" << current << "\n";
        }
        current = strtok(NULL, &m_delimiter);
    }
}

/** \brief Private function to process a Fixed track definition line in a wig file. Sets m_info details
 * \return
 */
void uParser::_processFixedWigLine(std::stringstream & curSStream)
{
    const std::string SPANSYMBOL="span=";
    const std::string STEPSYMBOL="step=";
    const std::string CHROMSYMBOL="chrom=";
    const std::string STARTSYMBOL="start=";

    auto pWigInfo =    std::dynamic_pointer_cast<wigInformation>(m_info);

    std::string chrom, strstart, strsspan;
    int curStart;
    /**< Chrom */
    curSStream >> chrom;
    if (chrom.find(CHROMSYMBOL)==std::string::npos)
        throw Parser_missing_mandatory_values()<<string_error("Missing chrom value in wig track definition");
    chrom=chrom.substr(chrom.find(CHROMSYMBOL)+CHROMSYMBOL.size());
    /**< Start */
    curSStream>>strstart;
    if (strstart.find(STARTSYMBOL)==std::string::npos)
        throw Parser_missing_mandatory_values()<<string_error("Missing start value in wig track definition");;

    strstart=strstart.substr(strstart.find(STARTSYMBOL)+STARTSYMBOL.size());
    curStart=stoi(strstart);
    /**< Step */
    std::string step;
    curSStream >>step;
    //Format definition is not sure if step is mandatory, so we will pretend it is not...
    int curStep=1;
    if (step.size())
    {
        if (step.find(STEPSYMBOL)==std::string::npos)
            throw  Parser_missing_mandatory_values()<<string_error("Missing step value in wig track definition");;;

        step=step.substr(step.find(STEPSYMBOL)+STEPSYMBOL.size());

        curStep=stoi(step);
    }
    /**< Optional Span parameter */
    int curSpan=1;
    std::string span;
    curSStream>>span;
    if(span.size())
    {
        /**< If not, fail again*/
        if (span.find(SPANSYMBOL)==std::string::npos)
            throw Parser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");

        span=span.substr(span.find(SPANSYMBOL)+SPANSYMBOL.size());
        curSpan=stoi(span);

    }
    //Nothing threw, modify values
    //Specification say's you cannot change spans
    if ((pWigInfo->getSpan()!=-1)&&(pWigInfo->getSpan()!=curSpan))
        throw Parser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file. Specification forbids changin span in dataset \n");
    pWigInfo->setStepType(wigInformation::stepType::FIXED);
    pWigInfo->setChrom(chrom);
    pWigInfo->setStep(curStep);
    pWigInfo->setSpan(curSpan);
    pWigInfo->setCurPos(curStart);
}

void uParser::_processVariabledWigLine(std::stringstream & curSStream)
{
    const std::string SPANSYMBOL="span=";
    const std::string CHROMSYMBOL="chrom=";


    auto pWigInfo =  std::dynamic_pointer_cast<wigInformation>(m_info);

    std::string chrom;
    std::string span;
    int curSpan=0;
    /**< Chrom tag */
    curSStream>> chrom;
    /**< If invalid chrom header */
    if (!(chrom.size())||((chrom.find(CHROMSYMBOL)==std::string::npos)))
        throw Parser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");
    chrom=chrom.substr(chrom.find(CHROMSYMBOL)+CHROMSYMBOL.size());
    /**< Optional Span parameter */
    curSpan=1;
    curSStream >> span;
    if(span.size())
    {
        /**< If good, yay, if not, fail again*/
        if (span.find(SPANSYMBOL)==std::string::npos)
            throw Parser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");;
        //   int testpos= span.find(SPANSYMBOL);
        span=span.substr(span.find(SPANSYMBOL)+SPANSYMBOL.size());
        curSpan=stoi(span);
    }
    if ((pWigInfo->getSpan()!=-1)&&(pWigInfo->getSpan()!=curSpan))
        throw Parser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file. Specification forbids changin span in dataset \n");
    pWigInfo->setStepType(wigInformation::stepType::VARIABLE);
    pWigInfo->setChrom(chrom);
    pWigInfo->setSpan(curSpan);
}

uToken uParser::_getNextEntryWig()
{
bool foundToken=false;
    try
    {
        while(foundToken==false){
            char line[4096];
            if (m_pIostream->getline(line, 4096))
            {
                std::stringstream ss;
                ss << line;

                std::string cur_token;
                ss>>cur_token;

                if((cur_token=="variableStep")||(cur_token=="fixedStep"))
                {
                    if (m_pIostream->eof())
                          throw end_of_file_throw()<<string_error("Badly formed wig file, track definition line with no entries following \n");
                    if (cur_token=="variableStep")
                    {
                        _processVariabledWigLine(ss);
                    }
                    else
                    {
                        _processFixedWigLine(ss);
                    }
                    //Process next line
                }
                else
                {

                    //If not eof
                    auto pWigInfo =  std::dynamic_pointer_cast<wigInformation>(m_info);
                    int end_pos;
                    int start_pos =0;
                    float score=0.0f;

                    switch (pWigInfo->getStepType())
                    {
                    case wigInformation::stepType::NA:

                        throw Parser_missing_mandatory_values()<<string_error("No declaraction line in Wig, error parsing \n");
                        break;
                    case  wigInformation::stepType::FIXED:
                        //curReg.setChr(curChr);
                        // curReg.setStartEnd(curStart,curStart+curSpan);
                        //curReg.setCount(utility::stringToInt(firstToken));
                        // curStart+=curStep;
                        break;
                    case   wigInformation::stepType::VARIABLE:
                        start_pos=stoi(cur_token);
                        ss >> cur_token;
                        score=stof(cur_token);
                        end_pos= start_pos+pWigInfo->getSpan();
                        break;
                    }
                    std::stringstream token_infos;

                    token_infos << "CHR\t" << pWigInfo->getChrom() << "\n";
                    token_infos << "START_POS\t" << start_pos << "\n";
                    token_infos << "END_POS\t" << end_pos << "\n";
                    token_infos << "SCORE\t" << score << "\n";
                    foundToken=true;
                    return uToken(token_infos);
                }
            }
            else
            {
    #ifdef DEBUG
                std::cerr << "Reached end of file." << std::endl;
    #endif
                end_of_file_throw e;
                e << string_error("Reached end of file.");
                throw e;
            }
        }
    }
    catch(invalid_uToken_throw& e)
    {
        throw e;
    }
    catch(Parser_missing_mandatory_values& e)
    {
        throw e;
    }
    catch(std::exception & e)
    {
        throw e;
    }

}

void uParser::_customParserValidateFields(const std::vector<std::string>& fieldsNames)
{
    /**< Must have at least 2 fields */
    if (fieldsNames.size() < 3)
    {
        customParser_missing_mandatory_values e;
        e << string_error("Custom file fields description is too short, must have at least 3 values: CHR and START_POS and a way to infer END_POS.\n");
        throw e;
    }

    /**< Check if mandatory fields are present */
    if (!_paramExists("CHR", fieldsNames))
    {
        customParser_missing_mandatory_values e;
        e << string_error("Mandatory field is missing: CHR\n");
        throw e;
    }
    if (!_paramExists("START_POS", fieldsNames))
    {
        customParser_missing_mandatory_values e;
        e << string_error("Mandatory field is missing: START_POS\n");
        throw e;
    }

    /**< We need to be able to infer END_POS either directly or indirectly (with a sequence or cigar score) */
    if (!_paramExists("END_POS", fieldsNames) && !_paramExists("SEQUENCE", fieldsNames) && !_paramExists("CIGAR", fieldsNames))
    {
        customParser_missing_mandatory_values e;
        e << string_error("We must be able to infer END_POS directly or indirectly (with SEQUENCE or CIGAR)\n");
        throw e;
    }
}

bool uParser::_paramExists(const std::string& name, const std::vector<std::string>& list) const
{
    std::vector<std::string>::const_iterator it;
    it = find(list.begin(), list.end(), name);
    if (it == list.end())
    {
        return false;
    }
    return true;
}

void uParser::_customParserCopyFields(const std::vector<std::string>& fieldsNames)
{
    /**< Only copy fields that have a matching token_param value, otherwise add NA */
    for(size_t i = 0; i < fieldsNames.size(); i++)
    {
        if (uToken::checkParam(fieldsNames[i]) == true)
        {
            m_customFieldNames.push_back(fieldsNames[i]);
        }
        else
        {
            m_customFieldNames.push_back("NA");
        }
    }
}
