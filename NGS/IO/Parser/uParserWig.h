// ***************************************************************************
// uParserWig.h (c) 2013
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



#ifndef UPARSERWIG_H_INCLUDED
#define UPARSERWIG_H_INCLUDED

#include "uParserBase.h"
#include "../../utility/utility.h"
#include <iostream>
//#include "uParserFactory.h"
namespace NGS
{

class uParserWig: public uParserBase
{
    class wigInformation
    {
    public:

        ~wigInformation() {};
        enum class stepType
        {
            FIXED,VARIABLE, NA
        };

        stepType getStepType()
        {
            return m_stepType;
        };
        void setStepType(stepType p_step)
        {
            m_stepType=p_step;
        };

        std::string getChrom()
        {
            return m_chrom;
        };
        void setChrom(std::string p_chrom)
        {
            m_chrom=p_chrom;
        };
        void setSpan(long int p_span)
        {
            m_span=p_span;
        };
        long int getSpan()
        {
            return m_span ;
        };
        void setStep(long int p_step)
        {
            m_step=p_step;
        };
        long int getStep()
        {
            return m_step ;
        };

        void setCurPos(long int p_curPos)
        {
            m_curPos=p_curPos;
        };
        long int getCurPos()
        {
            return m_curPos ;
        };
    private :
        int m_step=-1;
        long int m_curPos=0;
        stepType m_stepType=stepType::NA;
        std::string m_chrom="";
        long int m_span=-1;
    };

public :
    uParserWig();
    ~uParserWig();
    virtual void init(const std::string& filename, bool header = false);
    virtual void init(std::istream* stream, bool header = false);

    virtual uToken getNextEntry();

     static uParserBase * Create() { return new uParserWig(); }

private:
  //  static DerivedParserRegister<uParserWig> reg;
    wigInformation m_Info;
    void _processFixedWigLine(std::vector<std::string> & curSStream);
    void _processVariabledWigLine(std::vector<std::string> & curSStream);
	std::vector<std::string> m_tokens;
};

}
#endif // UPARSERWIG_H_INCLUDED
