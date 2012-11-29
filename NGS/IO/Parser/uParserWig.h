#ifndef UPARSERWIG_H_INCLUDED
#define UPARSERWIG_H_INCLUDED

#include "uParserBase.h"
#include <iostream>

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

private:
    static DerivedParserRegister<uParserWig> reg;
    wigInformation m_Info;
    void _processFixedWigLine(std::vector<std::string> & curSStream);
    void _processVariabledWigLine(std::vector<std::string> & curSStream);
	std::vector<std::string> m_tokens;
};

}
#endif // UPARSERWIG_H_INCLUDED
