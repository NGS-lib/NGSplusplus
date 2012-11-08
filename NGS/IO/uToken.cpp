#include "uToken.h"
#include "../utility/utility.h"
namespace NGS {

/** \brief uToken constructor.
 * \param istream& paramList: a stream containing all the parameters in the format: <token_param>\t<value>\n<token_param>\t<value>\n...
 * \exception invalid_value_throw when using a valid token param, but with an incorrect value (i.e.: a negative value for END_POS)
 * \exception invalid_uToken_throw when the token is in an invalid state, even if all the value are valid by themselves (i.e.: START_POS > END_POS)
 */
uToken::uToken(std::istream& paramList, bool customValues)
{
    m_customValues = customValues;
    char line[1024];
    bool custom = false;
    while (paramList.getline(line, 1024))
    {
        /**< Fetch name */
        token_param name;
        std::string value;
        std::string customName;
        try
        {
            std::stringstream ss;
            ss << line;
            ss >> name >> value;
            custom = false;
        }
        catch (invalid_token_param_throw& e)
        {
        /**< Only check for custom values is user wants them */
            if (m_customValues == true)
            {
                std::stringstream ss;
                ss << line;
                ss >> customName >> value;
                custom = true;
            }
        }

        /**< Fetch value */
        if (custom == false)
        {
            try
            {
                _setParam(name, value);
            }
            catch (invalid_value_throw& e)
            {
                throw e;
            }
        }
        else
        {
            if (m_customValues == true)
            {
                _setParamCustom(customName, value);
            }
        }
    }

    /**< Check if uToken is in a valid state  */
    try
    {
        _validateToken();
    }
    catch (invalid_uToken_throw& e)
    {
        addStringError(e,"Throwing in uToken constructor(std::istream paramList)");
        throw e;
    }
}

/** \brief Copy assignment operator
 */
uToken& uToken::operator=(uToken const& assign_from)
{
    if (this == &assign_from) return *this;
    m_params = assign_from.m_params;
    return *this;
}

/** \brief Get the value associated with a token_param.
 * \param token_param name: the name in token_param format.
 * \exception param_not_found when the param is not setted in the token.
 */
std::string uToken::getParam(token_param name) const
{
    if(isParamSet(name))
    {
        return (m_params.find(name)->second);
    }
    else
    {
        param_not_found e;
        std::string error = "Tried to getParam that is not set: " + _convertTokenParamToString(name) + "\n";
        addStringError(e,error);
        throw e;
    }
}

/** \brief Get the value associated with a parameter. Can be a token_param or a custom value.
 * \param const std::string& name: the name in string format.
 * \exception param_not_found when the param is not setted in the token.
 */
std::string uToken::getParam(const std::string& name) const
{
    try
    {
        token_param tokenName = _convertStringToTokenParam(name);
        if (isParamSet(tokenName))
        {
            return m_params.find(tokenName)->second;
        }
    }
    catch (invalid_token_param_throw& e)
    {
		if (m_customValues == false)
		{
			throw e;
		}
    }
        /**< If not a valid token_param, it can still be a custom value */
    if (isParamSet(name))
    {
        return m_customParams.find(name)->second;
    }
    else
    {
        param_not_found e;
        std::string error = "Tried to getParam that is not set: " + name + "\n";
        addStringError(e, error);
        throw e;
    }
}

/** \brief Set a param only if it's format is valid
* \param token_param& name: type of param
* \param const std::string& value: the value of the parameter
*/
void uToken::_setParam(const token_param& name, const std::string& value)
{
// TODO: What should we do in case the name is already setted? Overwrite silently, warning or error?
    try
    {
        if (_validateParam(name, value) == false)
        {
            invalid_value_throw e;
            e << string_error("Invalid token type is : "+this->_convertTokenParamToString(name)+"\n Invalid value is :"+value+"\n");
            throw e;
        }
        _postProcessParam(name, value);
        m_params[name] = value;
    }
    catch(invalid_value_throw &e)
    {
        addStringError(e,"Throwing in _setParam()");
        throw e;
    }
}

void uToken::_setParamCustom(const std::string& name, const std::string& value)
{
    if (m_customValues == true)
    {
        m_customParams[name] = value;
    }
}

/** \brief Check is the param is set.
 * \param const token_param& name: the name of the param in token_param format.
 * \return true if the param is set, otherwise return false.
 */
bool uToken::isParamSet(const token_param& name) const
{
    return m_params.count(name);
}

/** \brief Check is the param is set.
 * \param const std::string& name: the name of the param in string format.
 * \return true if the param is set, otherwise return false.
 */
bool uToken::isParamSet(const std::string& name) const
{
    try
    {
        token_param tokenName = _convertStringToTokenParam(name);
        return isParamSet(tokenName);
    }
    catch (invalid_token_param_throw& e)
    {
        if (m_customValues == true)
        {
            return m_customParams.count(name);
        }
    }
    return false;
}

/** \brief Operate any necessary post_processing on our values. Note this may include setting values to the token
 *
 * \param name token_param& Type of parameter
 * \param value const std::string& Value of parameter
 * \return void
 *
 */
void uToken::_postProcessParam(const token_param& name, const std::string& value) {
    switch(name) {
    case token_param::START_POS:
    case token_param::END_POS:
        break;
    case token_param::STRAND:
        break;
    case token_param::MAP_SCORE:
        break;
    case token_param::PHRED_SCORE:
        break;
    case token_param::CIGAR:
        _postProcCigar(value);
        break;
    case token_param::SEQUENCE:
        _postProcSequence(value);
        break;
    case token_param::FLAGS:
        _postProcFlag(value);
        break;
    case token_param::CHR:
        break;
    case token_param::SEQ_NAME:
        break;
    case token_param::SCORE:
        break;
    default:
        break;
    }
}

/** \brief Validate the parameter. Some tests are general, other are specific for each parameter.
 * \param token_param& name: type of param
 * \param const std::string& value: the value of the parameter
 */
bool uToken::_validateParam(const token_param& name, const std::string& value) const
{
    /**< In every case, an empty string is an invalid value */
    if (value.size() == 0)
    {
        return false;
    }
    //TODO numerical checks?
    /**< Then we do basic check for every type of parameters */
    switch(name)
    {
    case token_param::START_POS:
    case token_param::END_POS:
        return _posIsValid(value);
    case token_param::STRAND:
        return _strandIsValid(value);
    case token_param::MAP_SCORE:
        return _mapScoreIsValid(value);
    case token_param::PHRED_SCORE:
        return true;
    case token_param::CIGAR:
        return _cigarIsValid(value);
    case token_param::SEQUENCE:
        return _sequenceIsValid(value);
    case token_param::FLAGS:
        return _seqFlagsIsValid(value);
    case token_param::CHR:
        return true;
    case token_param::SEQ_NAME:
        return true;
    case token_param::SCORE:
         return _scoreIsValid(value);

    default:
        return false;
    }
}

void uToken::_validateToken()
{
    try {
        _checkMandatoryValues();
        _validateStartEnd();
        _validateSequenceCigar();
        _validateSequencePhred();
    }
    catch(invalid_uToken_throw &e) {
        throw e;
    }
}

void uToken::_checkMandatoryValues() const
{
    std::string str_chr = getParam(token_param::CHR);
    std::string str_start_pos = getParam(token_param::START_POS);
    std::string str_end_pos = getParam(token_param::END_POS);
    if (str_chr.size() == 0 || str_start_pos.size() == 0 || str_end_pos.size() == 0)
    {
        std::string error = "CHR, START_POS and END_POS are mandatory.\n";
        error += "CHR: " + str_chr + "\n";
        error += "START_POS: " + str_end_pos + "\n";
        error += "END_POS: " + str_start_pos + "\n";
        _throwInvalidToken(error);
    }
}




void uToken::_validateStartEnd() const
{
    std::string str_start_pos = getParam(token_param::START_POS);
    int int_start_pos = atoi(str_start_pos.c_str());
    std::string str_end_pos = getParam(token_param::END_POS);
    int int_end_pos = atoi(str_end_pos.c_str());
    if (int_start_pos > int_end_pos)
    {
        std::string error = "Invalid START_POS/END_POS values. \n";
        error += "START_POS:" + str_start_pos + ". END_POS: " + str_end_pos + ".";
        _throwInvalidToken(error);
    }
}

void uToken::_validateSequenceCigar() const {
    if (isParamSet(token_param::SEQUENCE))
    {
        std::string str_sequence;
        str_sequence = getParam(token_param::SEQUENCE);
        if (isParamSet(token_param::CIGAR))
    {
            std::string str_cigar = getParam(token_param::CIGAR);
            /**< Fetch all numerical values and sum them  */
            size_t cigar_size = 0;
            std::stringstream ss;
            for (size_t i = 0; i < str_cigar.size(); i++)
        {
                if (_isDigit(str_cigar[i]) == true)
        {
                    ss << str_cigar[i];
                }
                else {
                    int value;
                    switch(str_cigar[i])
            {
                    case 'M':
                    case 'I':
                    case 'S':
                    case '=':
                    case 'X':
                        ss >> value;
                        break;
                    default:
                        value = 0;
                        break;
                    }
                    cigar_size += value;
                    /**< Clear the stream */
                    ss.clear();
                    ss.str(std::string());
                }
            }
            /**< Check if sequence length and cigar sum match  */
            if (str_sequence.size() != cigar_size)
        {
                std::stringstream seq_size_stream;
                seq_size_stream << str_sequence.size();
                std::stringstream cigar_size_stream;
                cigar_size_stream << cigar_size;
                std::string error = "Sequence length and sum of cigar M/I/S/=/X values do not match.\n";
                error += "Sequence length: " + seq_size_stream.str();
                error += ". Sum of cigar values: " + cigar_size_stream.str() + "\n";
                error += "Cigar string is : "+str_cigar+"\n";
                _throwInvalidToken(error);
            }
        }
    }
}

void uToken::_validateSequencePhred() const
{
    if (isParamSet(token_param::SEQUENCE))
    {
        std::string str_sequence;
        str_sequence = getParam(token_param::SEQUENCE);
        if (isParamSet(token_param::PHRED_SCORE))
    {
            std::string str_phred = getParam(token_param::PHRED_SCORE);
            if (str_phred.size() != 0)
        {
                /**< If only 1 element and *, then not stored */
                if(! ((str_phred.size()==1) &&(str_phred.at(0)=='*') ))
        {
                    if (str_sequence.size() != str_phred.size())
            {
                        std::stringstream seq_size;
                        seq_size << str_sequence.size();
                        std::stringstream phred_size;
                        phred_size << str_phred.size();
                        std::string error = "Sequence and phred score length does not match.\n";
                        error += "Sequence size: " + seq_size.str() + ". Phred size: " + phred_size.str() + ".";
                        _throwInvalidToken(error);
                    }
                }
            }
        }
    }
}

void uToken::_throwInvalidToken(const std::string& baseErrorMessage) const
{
    invalid_uToken_throw e;
    e << string_error("Invalid token. " + baseErrorMessage);
    throw e;
}

/** \brief check if START_POS or END_POS is valid.
 * \param const std::string& value: the value of the pos entry.
 */
bool uToken::_posIsValid(const std::string& value) const
{
    std::stringstream ss;
    ss << value;
    /**< Check if the value is a long int greater than 0 */
    long int pos = 0;
    ss >> pos;
    if (pos <= 0)
    {
        return false;
    }
    /**< Finally, check if there is garbage at the end of value */
    return _isStreamEmpty(ss);
}

/** \brief Check if strand value is valid.
 * \param const std::string& value: the value of the strand entry.
 */
bool uToken::_strandIsValid(const std::string& value) const
{
    /**< Check if size is equal to 1 */
    if (value.size() != 1)
    {
        return false;
    }
    /**< Check if value is either '+' or '-' */
    if (value != "+" && value != "-")
    {
        return false;
    }
    return true;
}
/** \brief Check if map score value is valid (between 0 and 255 incl.).
 * \param const std::string& value: the value of the map score entry.
 */
bool uToken::_mapScoreIsValid(const std::string& value) const
{
    std::stringstream ss;
    ss << value;
    /**< Check if value is between */
    int map = -1;
    ss >> map;
    if (map < 0 || map > 255)
    {
        return false;
    }
    /**< Finally, check if there is garbage at the end of value */
    return _isStreamEmpty(ss);
}

/** \brief Check is sequence is valid.
 * \param const std::string& value: the value of the sequence entry.
 */
bool uToken::_sequenceIsValid(const std::string& value) const
{
    /**< Check if the first value is '*' and the size is 1 */
    if (value[0]== '*' && value.size() > 1)
    {
        return false;
    }
    /**< Otherwise, check if all the characters are in a valid format */
    for (size_t i = 0; i < value.size(); i++)
    {
        switch (value[i])
    {
        case 'a':
            break;
        case 'A':
            break;
        case 'c':
            break;
        case 'C':
            break;
        case 'g':
            break;
        case 'G':
            break;
        case 't':
            break;
        case 'T':
            break;
        case 'n':
            break;
        case 'N':
            break;
        case '.':
            break;
        case '=':
            break;
        case '*':
            break;
        default:
            return false;
        }
    }
    return true;
}

/** \brief Check if flag value is valid (between 0 and 65235 incl.).
 * \param const std::string& value: the value of the flag entry.
 */
bool uToken::_seqFlagsIsValid(const std::string& value) const
{
    std::stringstream ss;
    ss << value;
    /**< Check if value is between 0 and 65235 incl. */
    int flags = -1;
    ss >> flags;
    if (flags < 0 || flags > 65535)
    {
        return false;
    }
    /**< Finally, check if there is garbage at the end of value */
    return _isStreamEmpty(ss);
}
//TODO extra validation
/** \brief Check if cigar value is valid.
 * \param const std::string& value: the value of the cigar entry.
 */
bool uToken::_cigarIsValid(const std::string& value) const
{
    /**< Check if the first value is '*' and the size is 1 */
    if (value[0] == '*' && value.size() > 1)
    {
        return false;
    }
    /**< First item of value has to be a digit.  */
    if (_isDigit(value[0]) != true)
    {
        return false;
    }
    /**< Make sure there is not two alphabetical character in a row */
    bool lastValue = true;
    for (size_t i = 0; i < value.size(); i++)
    {
        bool currentValue = _isDigit(value[i]);
        if (lastValue == false && currentValue == false)
    {
            return false;
        }
        lastValue = currentValue;
    }
    /**< Must end with alphabetical value */
    if (_isDigit(value[value.size()-1]) == true)
    {
        return false;
    }
    return true;
}


bool uToken::_scoreIsValid(const std::string& value) const
{
    return true;
}


/** \brief Check if a value of a cigar score is a legal character.
 * \param char value: The value to test.
 */
bool uToken::_cigarValueIsValid(char value) const
{
    switch(value)
    {
    case 'M':
    case 'I':
    case 'D':
    case 'N':
    case 'S':
    case 'H':
    case 'P':
    case 'X':
    case '=': return true;
    default: return false;
    }
}

/** \brief Check if a istream is empty.
 * \param const std::istream& stream: the istream to check.
 */
bool uToken::_isStreamEmpty(const std::istream& stream) const
{
    if (stream.rdbuf()->in_avail() != 0)
    {
        return false;
    }
    return true;
}

/**< Post processing methods */

/** \brief Sequence length can be used to infer END_POS
 * \param sequence: The nucleotidic sequence of the entry.
 */

void uToken::_postProcSequence(const std::string& sequence) {
    if (!(isParamSet(token_param::END_POS))) {
        auto start_pos = std::stoi(getParam(token_param::START_POS));
        _setParam(token_param::END_POS, std::to_string(sequence.size()+start_pos-1));
    }
}

/** \brief Sam Flag can contain multiple settings
 *
 * \param value
 * \return void
 *
 */
void uToken::_postProcFlag(const std::string& flag)
{

}
void uToken::_postProcCigar(const std::string& cig) {
    /**< If END_POS is not set, calculate it's value from cigar score and set it */
    if (!(isParamSet(token_param::END_POS))) {
        try {
            int curPos=0;
            std::string substr;
            int size=0;
            for (unsigned int i=0; i< cig.size(); i++) {
                /**< If isAlpha then check previous numbers */
                if (isalpha(cig.at(i))) {
                    /**< If a count value */
                    char temp;
                    temp = cig.at(i);
                    if ((temp=='M')||(temp=='I')||(temp=='S')||(temp=='X')||(temp=='+')) {
                        substr= cig.substr(curPos, (i-curPos));
                        size+= atoi(substr.c_str());
                    }
                    curPos=(i+1);
                }
            }
            auto start_pos=std::stoi(getParam(token_param::START_POS));
            _setParam(token_param::END_POS, std::to_string(start_pos+(size-1) ));
        }
        catch(uToken_exception_base &e) {
            addStringError(e, "Throwing, in _postProcCigar, unable to set END_POS as START_POS not set");
            throw e;
        }
    }
}


std::string uToken::_convertTokenParamToString(const token_param& token) const
{
    std::stringstream ss;
    ss << token;
    return ss.str();
}

token_param uToken::_convertStringToTokenParam(const std::string& name) const
{
    std::stringstream ss;
    ss << name;
    token_param tokenParam;
    ss >> tokenParam;
    return tokenParam;
}

} // End of namespace NGS