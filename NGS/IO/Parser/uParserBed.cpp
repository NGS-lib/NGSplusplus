#include "uParserBed.h"

namespace NGS {

/** \brief Default constructor (not used directly).
 */
uParserBed::uParserBed(): uParserBase()
{
}

/** \brief Destructor.
 */
uParserBed::~uParserBed()
{
}

/** \brief Initialize the uParserBed object (open file and parse header).
 * \param const std::string& filename: name of the bed file to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserBed::init(const std::string& filename, bool header)
{
    uParserBase::init(filename, header);
    if (header == true)
    {
        _parseHeader();
    }
    m_headerParsed = true;
    m_delimiter = '\t';
}

/** \brief Initialize the uParserBed object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserBed::init(std::istream* stream, bool header)
{
    uParserBase::init(stream, header);
    m_headerParsed = true;
    m_delimiter = '\t';
}

/** \brief Produce a token with next entry in the file/stream.
 * \return uToken containing the infos of the next entry.
 */
uToken uParserBed::getNextEntry()
{
    char line[4096];

    if (m_pIostream->getline(line, 4096))
    {
        m_rawString=line;
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
    std::cerr <<"Fatal error in getNextEntry() from Bed parser, should not reach here." <<std::endl;
    abort();
}

void uParserBed::_convertLineToTokenInfosBed(char* line, std::stringstream& token_infos)
{
    char line_cpy[4096];
    strcpy(line_cpy, line);
    std::stringstream ss;
    ss << line_cpy;
    if (m_numberOfColumn == 0 && m_headerParsed == true)
    {
        int numberOfColumn = _countColumns(line);
        _validateColumnNumber(numberOfColumn);
        m_numberOfColumn = numberOfColumn;
    }
    std::string chr = _getNextEntry(line_cpy);
    std::string start_pos = _getNextEntry(line_cpy);
    std::string end_pos = _getNextEntry(line_cpy);
    std::string score;
    std::string strand;
    std::string extra;
    /**< BED3 to BED6 are supported */
    token_infos << "CHR\t" << chr << "\n";
    token_infos << "START_POS\t" << start_pos << "\n";
    token_infos << "END_POS\t" << end_pos << "\n";
    if (m_numberOfColumn > 3) {
        std::string seq_name = _getNextEntry(line_cpy);
        token_infos << "SEQ_NAME\t" << seq_name << "\n";
    }
    if (m_numberOfColumn > 4) {
        score = _getNextEntry(line_cpy);

        if (score!="."){
        token_infos << "SCORE\t" << score << "\n"; }

	}
    if (m_numberOfColumn > 5)
    {
        strand = _getNextEntry(line_cpy);
        if (strand!=".")
        token_infos << "STRAND\t" << strand << "\n";
    }
    if (m_numberOfColumn > 6) {
	extra = line_cpy;
	token_infos << "EXTRA\t" << extra << "\n";
    }
}

int uParserBed::_countColumns(char* line) const
{
    int count = 1;
    for (int i = 0; line[i] != '\0' && line[i] != '\n'; i++)
    {
        if (line[i] == '\t')
        {
            count++;
        }
    }
    return count;
}

void uParserBed::_validateColumnNumber(int numberOfColumn) const
{
    if (numberOfColumn < 3 || numberOfColumn > 6)
    {
        uParserBed_invalid_number_of_columns e;
        e << string_error("uParserBed: Invalid number of columns.");
        throw e;
    }
}

void uParserBed::_parseHeader()
{
    std::string unformatedHeader;
    /**< We parse a line at a time until we get a valid token, then we return the token back in the stream */
    do
    {
        char line[4096];
        if (m_pIostream->getline(line, 4096))
        {
            std::stringstream token_infos;
            _convertLineToTokenInfosBed(line, token_infos);
            try
            {
                uToken token(token_infos);
                /**< When we have a valid token, we consider that all the header have beed fetched */
                /**< We put token back into stream, and exit function */
                _pushBackLine(line);
                m_headerParsed = true;
            }
            /**< If there is a throw, we are not at the first valid token line, so we add line to header object */
            catch(invalid_uToken_throw& e) { }
            catch(invalid_value_throw& e) { }
            if (m_headerParsed == false)
            {
                std::string s(line);
                unformatedHeader += s;
                unformatedHeader += '\n';
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
    } while(m_headerParsed == false);
    m_headerData.setUnformatedHeader(unformatedHeader);
}

void uParserBed::_pushBackLine(char* line)
{
    m_pIostream->putback('\n');
    std::string str_line(line);
    for (int i = (int)(str_line.size()) - 1; i >= 0; i--)
    {
        m_pIostream->putback(line[i]);
    }
}

std::string uParserBed::_getNextEntry(char* line)
{
    char toReturnChar[4096];
    int i = 0;

    /**< Fetch the next entry */
    do
    {
        toReturnChar[i] = line[i];
        i++;
    } while (line[i] != m_delimiter && line[i] != '\0');
    toReturnChar[i] = '\0';

    /**< Remove fetched entry from line */
    int j = 0;
    if (line[i] != '\0')
    {
        i++;
        do
        {
            line[j] = line[i];
            i++;
            j++;
        } while (line[i] != '\n' && line[i] != '\0');
    }
    line[j] = '\0';

    std::string toReturn(toReturnChar);
    return toReturn;
}

//DerivedParserRegister<uParserBed> uParserBed::reg("BED");
} // End of namespace NGS
