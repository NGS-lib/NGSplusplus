// ***************************************************************************
// uHeader.h (c) 2013
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


#ifndef UHEADER_H_INCLUDED
#define UHEADER_H_INCLUDED

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include "../uGeneException.h"
#include "../utility/utility.h"

namespace NGS {
/**< List of param is hard coded as strongly typed enum for extra safety. */
/**< This list has to be updated for every new param */
enum class header_param { CHR, CHR_SIZE,STEP_SIZE, SORT_TYPE};

/**< uHeader class, to keep track of information in header, formated or not */
/**< This is the class that takes care of data validation */
/**< All the data is saved in a map in string format */
class uHeader {
public:
    /** \brief Default constructor.
      */
    uHeader() ;

    /** \brief Add string to header as is.
      * \param const std::string& unformatedHeader: the string to copy to header
      */
    void setUnformatedHeader(const std::string& unformatedHeader) { m_unformatedHeader = unformatedHeader; };
    void _setParam(const header_param& name, const std::string& value);
    void _addToParam(const header_param& name, const std::string& value);

    /** \brief Fetch a param. Throw param_not_found if the param does not exist.
     * \param header_param& name: the name of the param we wish to get.
     */
    std::vector<std::string> getParamVector(header_param name) const;
    std::string getParam(header_param name) const;

    /** \brief Return the exact same header that was in the input data.
      */
    std::string getUnformatedHeader() const { return m_unformatedHeader;};

    /** \brief Check if there is a value associated with a given param.
          * \param header_param& name: name of the param to check.
      */
    bool isParamSet(const header_param& name) const { return m_params.count(name); };
    uHeader& operator=(uHeader const& assign_from);

    /** \brief Check if a string has a corresponding header_param value
     * \param const std::string& param: The param to check.
     * \return True if there is a matching header_param value, otherwise false.
     */
    static inline bool checkParam(const std::string& param) {

        return (param == "CHR"
             || param == "CHR_SIZE"
             || param == "STEPS_SIZE"
             || param == "SORT_TYPE"
              );
    }

private:
    std::map<header_param, std::vector<std::string>> m_params={};
    std::string m_unformatedHeader="";

    void _postProcessParam(const header_param& name, std::string& value);
    bool _validateParam(header_param name, const std::string& value);

    std::string _convertHeaderParamToString(const header_param& header) const;
    std::map<header_param, std::function<bool(const uHeader*,const std::string&)> > validate_func_map;
    std::map<header_param, std::function<void(uHeader*,std::string)>> post_func_map;


    bool _noValidate(const std::string& value)const {return true;};
    void _noPost(std::string) {};

    /**< Chromosome header parameters */
    bool _validateChrSize(const std::string& sizeString)const;
    bool _validateChrList(const std::string& chrString)const;
    /**< Wig */
    bool _valideStepSize(const std::string& sizeString)const;
}; // End of class Header




/**< Overloading of stream operator for header_param */
inline std::ostream & operator<<(std::ostream& Str, header_param name) {
    switch (name) {
    case header_param::CHR : return Str << "CHR";
    case header_param::CHR_SIZE : return Str << "CHR_SIZE";
    case header_param::STEP_SIZE : return Str << "STEP_SIZE";
    case header_param::SORT_TYPE : return Str << "SORT_TYPE";
    default: return Str << (int) name;
    }
}

/**< Overloading of stream operator for header_param */
inline std::string& operator<<(std::string& Str, header_param name) {
    switch (name) {
    case header_param::CHR : return Str+="CHR";
    case header_param::STEP_SIZE : return Str+="STEP_SIZE";
    case header_param::CHR_SIZE : return Str+="CHR_SIZE";
    case header_param::SORT_TYPE : return Str+="SORT_TYPE";
    default: return Str;
    }
}

inline std::istream& operator>>(std::istream &is, header_param& name) {
    std::string header;
    is >> header;
    if (header == "CHR") name = header_param::CHR;
    else if (header == "CHR_SIZE") name = header_param::CHR_SIZE;
    else if (header == "STEP_SIZE") name = header_param::STEP_SIZE;
    else if (header == "SORT_TYPE") name = header_param::SORT_TYPE;
    else {
        invalid_header_param_throw e;
        e << string_error(header);
        throw e;
    }
    return is;
}
} // End of namespace NGS
#endif // UHEADER_H_INCLUDED
