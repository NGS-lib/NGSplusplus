// ***************************************************************************
// utility.cpp (c) 2013
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


#include "utility.h"

//Modified from
//http://stackoverflow.com/questions/53849/how-do-i-tokenize-a-string-in-c
namespace utility{

    const std::string Tokenizer::DELIMITERS(" \t\n\r");
    Tokenizer::Tokenizer(const std::string& s) :
        m_offset(0),
        m_string(s),
        m_delimiters(DELIMITERS) {}

    Tokenizer::Tokenizer(const std::string& s, const std::string& delimiters) :
         m_offset(0),
        m_string(s),
        m_delimiters(delimiters) {}

    bool Tokenizer::NextToken()
    {
        return NextToken(m_delimiters);
    }

    bool Tokenizer::NextToken(const std::string& delimiters)
    {
        size_t i = m_string.find_first_not_of(delimiters, m_offset);
        if (std::string::npos == i)
        {
            m_offset = m_string.length();
            return false;
        }

        size_t j = m_string.find_first_of(delimiters, i);
        if (std::string::npos == j)
        {
            m_token = m_string.substr(i);
            m_offset = m_string.length();
            return true;
        }

        m_token = m_string.substr(i, j - i);
        m_offset = j;
        return true;
    }
     const std::string Tokenizer::GetToken() const
    {

        return m_token;

    }
}
