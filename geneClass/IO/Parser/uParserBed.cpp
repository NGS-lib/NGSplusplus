#include "uParserBed.h"

namespace NGS {

uParserBed::uParserBed(): uParserBase() 
{
}

uParserBed::~uParserBed() 
{
	if (m_dynamicStream == true) 
	{
		delete m_pIostream;
	}
	m_pIostream = NULL;
}

void uParserBed::init(const std::string& filename, bool header) 
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
		m_dynamicStream = true;
	}
	if (header == true)
	{
		_parseHeader();
	}
	m_headerParsed = true;
}

void uParserBed::init(std::iostream* stream, bool header) 
{
	m_pIostream = stream;
	m_dynamicStream = false;
	m_headerParsed = true;
}

void uParserBed::init(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header, char delimiter)
{
    throw ugene_exception_base()<<string_error("Invalid constructor call for Bed Format");
}

void uParserBed::init(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header, char delimiter)
{
    throw ugene_exception_base()<<string_error("Invalid constructor call for Bed Format");
}

uToken uParserBed::getNextEntry() 
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
	std::cerr <<"Fatal error in getNextEntry() from Bed parser, should not reach here." <<std::endl;
	abort();
}

void uParserBed::_convertLineToTokenInfosBed(char* line, std::stringstream& token_infos)
{
	std::stringstream ss;
	ss << line;
	if (m_numberOfColumn == 0 && m_headerParsed == true) 
	{
		m_numberOfColumn = _countColumns(line);
		_validateColumnNumber();
	}
	std::string chr;
	std::string start_pos;
	std::string end_pos;
	std::string score;
	std::string seq_name;
	std::string strand;
	ss >> chr >> start_pos >> end_pos >> seq_name; // >> score >> strand;
	token_infos << "CHR\t" << chr << "\n";
	token_infos << "START_POS\t" << start_pos << "\n";
	token_infos << "END_POS\t" << end_pos << "\n";
	token_infos << "SEQ_NAME\t" << seq_name << "\n";
	/**< If there is no score info, we have a Bed4. We don't add SCORE and STRAND to uToken. */
//	if (score.size() != 0)
	if (m_numberOfColumn == 6)
	{
		ss >> score >> strand;
		token_infos << "SCORE\t" << score << "\n";
		token_infos << "STRAND\t" << strand << "\n";
//		if (score.size() == 0 || strand.size() == 0)
//		{
//			uParser_missing_mandatory_values e;
////			Parser_missing_mandatory_values e;
//			e << string_error("SCORE and STRAND values are mandatory in BED6.");
//			throw e;
//		}
	}
}

int uParserBed::_countColumns(char* line) const
{
	int count = 1;
	for (int i = 0; line[i] != '\0' || line[i] != '\n'; i++)
	{
		if (line[i] == '\t')
		{
			count++;
		}
	}
	return count;
}

void uParserBed::_validateColumnNumber() const
{
	if (m_numberOfColumn != 4 && m_numberOfColumn != 6) 
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

DerivedParserRegister<uParserBed> uParserBed::reg("BED");
} // End of namespace NGS
