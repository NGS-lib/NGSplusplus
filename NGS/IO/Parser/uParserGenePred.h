// ***************************************************************************
// uParserGenePred.h (c) 2013
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


#ifndef _H_INCLUDED
#define _H_INCLUDED

#include "uParserBase.h"
namespace NGS
{

    class uParserGenePred : public uParserBase
    {
    public :
        uParserGenePred();
        ~uParserGenePred();

        virtual void init(const std::string& filename, bool header = false);
        virtual void init(std::istream* stream, bool header = false);
        uToken getNextEntry();
        static uParserBase * Create() { return new uParserGenePred(); }
    private:
        uToken _getTokenFromGenePredString(const std::string & line);
        void  _parseHeader();
        boost::xpressive::sregex UCSCGenePredRegex;
        boost::xpressive::sregex UCSCGenePredRegexInit;
       // const std::string GenePredregStringPriortest="(.*)";
        const std::string GenePredregStringPrior="^([.\\w_-]+)\t([.\\w_-]+)\t(\\+|\\-)\t(\\d+)\t(\\d+)\t(\\d+)\t(\\d+)\t(\\d+)\t(.*)";
        const std::string GenePredregStringPart1="^([.\\w_-]+)\t([.\\w_-]+)\t(\\+|\\-)\t(\\d+)\t(\\d+)\t(\\d+)\t(\\d+)\t(\\d+)";
        const std::string GenePredRegExon ="(\\d+),";

    };

} // End of namespace NGS

#endif // _H_INCLUDED
