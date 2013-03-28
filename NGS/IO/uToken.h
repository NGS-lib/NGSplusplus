// ***************************************************************************
// uToken.h (c) 2013
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
// ***************************************************************************

#ifndef UTOKEN_H_INCLUDED
#define UTOKEN_H_INCLUDED

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include "../uGeneException.h"
#include "../utility/utility.h"

namespace NGS
{
/**< List of param is hard coded as strongly typed enum for extra safety. */
/**< This list has to be updated for every new param */
enum class token_param
{
    CHR, START_POS, END_POS, STRAND, MAP_SCORE, PHRED_SCORE, CIGAR, SEQUENCE, SEQ_NAME,FLAGS ,SCORE,DENSITY,FEATURE_TYPE,SOURCE,PHASE,EXTRA,TEMPLATE_LENGTH,
    FEATURE_CLASS, FEATURE_ID,GROUP_TRANSCRIPT,GROUP_ID,OFFSET,
    CUST_1,CUST_2,CUST_3,CUST_4,CUST_5,CUST_6,CUST_7,CUST_8,CUST_9
};

enum class validate_type
{
    VALIDATE, NO_VALIDATE
};

/**< uToken class, to bridge parser and the library's class */
/**< This is the class that takes care of data validation */
/**< All the data is saved in a map in string format */
class uToken
{
public:
//    uToken(std::istream& paramList, bool customValues = false, validate_type validate=validate_type::VALIDATE );
    uToken(std::istream& paramList, validate_type validate=validate_type::VALIDATE);
    uToken(std::istream& paramList, bool customValues, validate_type validate=validate_type::VALIDATE);
    uToken();

    /** \brief Fetch a param. Throw param_not_found if the param does not exist.
     * \param token_param& name: the name of the param we wish to get.
     */
    std::string getParam(token_param name,unsigned int paramPos=0) const;
    std::string getParam(const std::string& name,unsigned int paramPos=0) const;

    void print(std::ostream &pOut)const;
    bool isParamSet(const token_param& name, int pos=0) const;
    int  paramCount(const token_param& name) const;
    bool isParamSet(const std::string& name, int pos=0) const;
    int  paramCount(const std::string& name) const;
    std::vector<std::string> getParamVector(token_param name) const;
    uToken& operator=(uToken const& assign_from);

    /** \brief Check if a string has a corresponding token_param value
     * \param const std::string& param: The param to check.
     * \return True if there is a matching token_param value, otherwise false.
     */
    static inline bool checkParam(const std::string& param)
    {
        return (param == "CHR"
                || param == "START_POS"
                || param == "END_POS"
                || param == "STRAND"
                || param == "MAP_SCORE"
                || param == "PHRED_SCORE"
                || param == "CIGAR"
                || param == "SEQUENCE"
                || param == "SEQ_NAME"
                || param == "FLAGS"
                || param == "DENSITY"
                || param == "SCORE"
                || param == "FEATURE_TYPE"
                || param == "FEATURE_CLASS"
                || param == "FEATURE_ID"
                || param == "FEATURE_TRANSCRIPT"
                || param == "SOURCE"
                || param == "PHASE"
                || param == "EXTRA"
		|| param == "TEMPLATE_LENGTH"
		|| param == "GROUND_TRANSCRIPT"
		|| param == "GROUND_ID"
		|| param == "OFFSET"
                );
    }

private:
    void _setParamNoValidate(const token_param& name, const std::string& value);
    std::map<token_param, std::vector<std::string>> m_params= {};
    bool m_customValues = false;
    std::map<std::string, std::string> m_customParams= {};


    void _initialize(std::istream& paramList, bool customValues, validate_type validate);

    void _setParam(const token_param& name, const std::string& value);
    void _setParamCustom(const std::string& name, const std::string& value);

    void _postProcessParam(const token_param& name, const std::string& value);
    bool _validateParam(const token_param& name, const std::string& value) const;
    void _validateToken();
    void _checkMandatoryValues() const;
    void _validateStartEnd() const;
    void _validateSequenceCigar() const;
    void _validateSequencePhred() const;
    void _throwInvalidToken(const std::string& baseErrorMessage) const;

    /**< Single parameter validation functions */
//    bool _chrIsValid(const std::string& name) const;
    bool _posIsValid(const std::string& name) const;
    bool _strandIsValid(const std::string& name) const;
    bool _mapScoreIsValid(const std::string& name) const;
//    bool _phredScoreIsValid(const std::string& name) const;
    bool _sequenceIsValid(const std::string& name) const;
//    bool _seqNameIsValid(const std::string& name) const;
    bool _scoreIsValid(const std::string& name) const;
    bool _seqFlagsIsValid(const std::string& name) const;
    bool _cigarIsValid(const std::string& name) const;

    bool _isDigit(char value) const
    {
        return std::isdigit(value);
    }
    //(value >= '0' && value <= '9'); }
    bool _cigarValueIsValid(char value) const;
    bool _isStreamEmpty(const std::istream& stream) const;

    void _postProcSequence(const std::string& sequence);
    void _postProcFlag(const std::string& flag);
    void _postProcCigar(const std::string& cig);

    std::string _convertTokenParamToString(const token_param& token) const;
    token_param _convertStringToTokenParam(const std::string& name) const;


    /**< OUr support formats have access to the Token as they validate themselves */

	friend class uParserBase;
	friend class uParserSam;
	friend class uParserWig;
	friend class uParserGFF;
	friend class uParserBedGraph;
	friend class uParserGTF;
	friend class uParserGFF2;
	friend class uParserGenePred;
	friend class uParserBAM;
}; // End of class Token

/**< Overloading of stream operator for token_param */
inline std::ostream & operator<<(std::ostream& Str, token_param name)
{
    switch (name)
    {
    case token_param::CHR:
        return Str << "CHR";
    case token_param::START_POS:
        return Str << "START_POS";
    case token_param::END_POS:
        return Str << "END_POS";
    case token_param::STRAND:
        return Str << "STRAND";
    case token_param::MAP_SCORE:
        return Str << "MAP_SCORE";
    case token_param::PHRED_SCORE:
        return Str << "PHRED_SCORE";
    case token_param::CIGAR:
        return Str << "CIGAR";
    case token_param::SEQUENCE:
        return Str << "SEQUENCE";
    case token_param::SEQ_NAME:
        return Str << "SEQ_NAME";
    case token_param::FLAGS:
        return Str << "FLAGS";
    case token_param::SCORE:
        return Str <<"SCORE";
    case token_param::DENSITY:
        return Str <<"DENSITY";
    case token_param::FEATURE_TYPE:
        return Str <<"FEATURE_TYPE";
    case token_param::FEATURE_CLASS:
        return Str <<"FEATURE_CLASS";
    case token_param::FEATURE_ID:
        return Str <<"FEATURE_ID";
    case token_param::GROUP_TRANSCRIPT:
        return Str <<"GROUP_TRANSCRIPT";
    case token_param::GROUP_ID:
        return Str <<"GROUP_ID";
     case token_param::PHASE:
        return Str <<"PHASE";
     case token_param::SOURCE:
        return Str <<"SOURCE";
     case token_param::EXTRA:
        return Str <<"EXTRA";
	case token_param::TEMPLATE_LENGTH:
        return Str <<"TEMPLATE_LENGTH";
    case token_param::OFFSET:
        return Str <<"OFFSET";

    default:
        return Str << (int) name;
    }
}

inline std::istream& operator>>(std::istream &is, token_param& name)
{
    std::string token;
    is >> token;
    if (token == "CHR") name = token_param::CHR;
    else if (token == "START_POS") name = token_param::START_POS;
    else if (token == "END_POS") name = token_param::END_POS;
    else if (token == "STRAND") name = token_param::STRAND;
    else if (token == "MAP_SCORE") name = token_param::MAP_SCORE;
    else if (token == "PHRED_SCORE") name = token_param::PHRED_SCORE;
    else if (token == "CIGAR") name = token_param::CIGAR;
    else if (token == "SEQUENCE") name = token_param::SEQUENCE;
    else if (token == "SEQ_NAME") name = token_param::SEQ_NAME;
    else if (token == "FLAGS") name = token_param::FLAGS;
    else if (token == "SCORE") name = token_param::SCORE;
    else if (token == "DENSITY") name = token_param::DENSITY;
    else if (token == "FEATURE_TYPE") name = token_param::FEATURE_TYPE;
    else if (token == "PHASE") name = token_param::PHASE;
    else if (token == "SOURCE") name = token_param::SOURCE;
    else if (token == "EXTRA") name = token_param::EXTRA;
    else if (token == "TEMPLATE_LENGTH") name = token_param::TEMPLATE_LENGTH;
    else if (token == "FEATURE_CLASS") name = token_param::FEATURE_CLASS;
    else if (token == "FEATURE_ID") name = token_param::FEATURE_ID;
    else if (token == "GROUP_ID") name = token_param::GROUP_ID;
    else if (token == "GROUP_TRANSCRIPT") name = token_param::GROUP_TRANSCRIPT;
    else if (token == "OFFSET") name = token_param::OFFSET;
    else
    {
        invalid_token_param_throw e;
        e << string_error("Failling in Token >> operator, missing param: "+token+"\n");
        throw e;
    }
    return is;
}
} // End of namespace NGS
#endif // UTOKEN_H_INCLUDED
