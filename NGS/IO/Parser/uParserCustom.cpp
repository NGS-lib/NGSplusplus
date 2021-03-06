// ***************************************************************************
// uParserCustom.cpp (c) 2013
// Alexei Nordell-Markovits : Sherbrooke University
// Charles Joly Beauparlant : Laval University
//
//       This file is part of the NGS++ library.
//
//    The NGS++ library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program (lgpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// **************************************************************************



#include "uParserCustom.h"

namespace NGS {

/** \brief Default constructor (not to be used directly, initialize through uParser instead)
 */
uParserCustom::uParserCustom(): uParserBed()
{
}

/** \brief Destructor.
 */
uParserCustom::~uParserCustom()
{
}

/** \brief Initialize the uParserCustom object (default initialization of uParserCustom is not valid).
 */
void uParserCustom::init(const std::string& filename, bool header)
{
    throw uParser_exception_base()<<string_error("Invalid constructor call for Custom Format");
}

/** \brief Initialize the uParserCustom object (default initialization of uParserCustom is not valid).
 */
void uParserCustom::init(std::istream* stream, bool header)
{
    throw uParser_exception_base()<<string_error("Invalid constructor call for Custom Format");
}

/** \brief Initialize the uParserCustom object.
 * \param const std::string& filename: the file to parse.
 * \param const std::vector<std::string>& fieldNames: the identification of columns in the file.
 * \param char delimiter: the delimitor between each field in a row (default: tabulation).
 */
void uParserCustom::init(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter)
{
    uParserBase::init(filename);
    m_delimiter = delimiter;
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
}

/** \brief Initialize the uParserCustom object.
 * \param std::istream* stream: the stream to parse.
 * \param const std::vector<std::string>& fieldNames: the identification of columns in the file.
 * \param char delimiter: the delimitor between each field in a row (default: tabulation).
 */
void uParserCustom::init(std::istream* stream, const std::vector<std::string>& fieldsNames, char delimiter)
{
    uParserBase::init(stream);
    m_delimiter = delimiter;
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
}

/** \brief Produce a token with next entry in the file/stream.
 * \return uToken containing the infos of the next entry.
 */
uToken uParserCustom::getNextEntry()
{
    char line[4096];
    if (m_pIostream->getline(line, 4096))
    {
        /**< We start by fetching the infos from the line */
        m_rawString=line;
        std::stringstream token_infos;

        _convertLineToTokenInfosCustom(line, token_infos);
        /**< We try to create a token with the infos that were fetched from the line */
        /**< If it doesn't work andit's the first token, we don't throw an error. Instead, we try again with another line */
        try
        {
            uToken token(token_infos, true);
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
    std::cerr <<"Fatal error in getNextEntry() from Custom parser, should not reach here." <<std::endl;
    abort();
}

void uParserCustom::_convertLineToTokenInfosCustom(char* line, std::stringstream& token_infos)
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
        if (m_customFieldNames[i]!="JUNK"){
            token_infos << m_customFieldNames[i] << "\t" << current << "\n";
        }
            current = strtok(NULL, &m_delimiter);

    }
}

void uParserCustom::_customParserValidateFields(const std::vector<std::string>& fieldsNames) const
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

bool uParserCustom::_paramExists(const std::string& name, const std::vector<std::string>& list) const
{
    std::vector<std::string>::const_iterator it;
    it = find(list.begin(), list.end(), name);
    if (it == list.end())
    {
        return false;
    }
    return true;
}

void uParserCustom::_customParserCopyFields(const std::vector<std::string>& fieldsNames)
{
    /**< Only copy fields that have a matching token_param value, otherwise add NA */
    for(size_t i = 0; i < fieldsNames.size(); i++)
    {
        m_customFieldNames.push_back(fieldsNames[i]);
    }
}

//DerivedParserRegister<uParserCustom> uParserCustom::reg("CUSTOM");
} // End of namespace NGS
