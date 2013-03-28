// ***************************************************************************
// uWriterCustom.cpp (c) 2013
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


#include "uWriterCustom.h"

namespace NGS
{


/** \brief Constructor to set default values (the user have to create the class through the writer class)
  */
uWriterCustom::uWriterCustom()
{
    m_delimiter = '\t';
    m_fieldsNames.push_back("CHR");
    m_fieldsNames.push_back("START_POS");
    m_fieldsNames.push_back("END_POS");
    m_fieldsNames.push_back("SEQ_NAME");

}

/** \brief Print the values of a token in Custom format in current file
  * \param const uToken& token: the token to print.
  */
void uWriterCustom::writeToken(const uToken& token)
{
    std::vector<std::string> values;
    for (size_t i = 0; i < m_fieldsNames.size(); i++)
    {
        if (uToken::checkParam(m_fieldsNames[i]) == true)
	{
            std::stringstream ss;
            ss << m_fieldsNames[i];
            token_param param;
            try
	    {
                ss >> param;
                if (token.isParamSet(param))
		{
                    values.push_back(token.getParam(param));
                }
                else
		{
                    values.push_back(".");
                }
            }
            catch (invalid_token_param_throw& e)
	    {
                values.push_back(".");
            }
        }
        else
	{
            values.push_back(".");
        }
        *m_pOstream << values[i];
        if (i != m_fieldsNames.size() - 1)
	{
            *m_pOstream << m_delimiter;
        }
    }
    *m_pOstream << std::endl;
}

//DerivedRegister<uWriterCustom> uWriterCustom::reg("CUSTOM");

} // End of namespace NGS
