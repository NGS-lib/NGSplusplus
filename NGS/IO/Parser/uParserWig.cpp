#include "uParserWig.h"

namespace NGS
{

uParserWig::uParserWig(): uParserBase()
{
}

uParserWig::~uParserWig(){

}

void uParserWig::init(const std::string& filename, bool header )
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
         char line[4096];
    if (m_pIostream->getline(line, 4096))
    {
        std::stringstream ss;
        ss << line;

        std::string cur_token;
        ss>>cur_token;

        if((cur_token=="variableStep")||(cur_token=="fixedStep"))
        {
            if (m_pIostream->eof())
                throw end_of_file_throw()<<string_error("Badly formed wig file, track definition line with no entries following \n");
           try {
            if (cur_token=="variableStep")
            {
                _processVariabledWigLine(ss);
            }
            else
            {
                _processFixedWigLine(ss);
            }
            }
            catch(uParser_missing_mandatory_values &e){
                throw e;
            }
        }
        else
            throw uParser_missing_mandatory_header()<<string_error("Missing header in wiggle file \n");

        }
    }

}
void uParserWig::init(std::iostream* stream, bool header )
{
    m_pIostream = stream;
    m_dynamicStream = false;
      char line[4096];
    if (m_pIostream->getline(line, 4096))
    {
        std::stringstream ss;
        ss << line;

        std::string cur_token;
        ss>>cur_token;

        if((cur_token=="variableStep")||(cur_token=="fixedStep"))
        {
            if (m_pIostream->eof())
                throw end_of_file_throw()<<string_error("Badly formed wig file, track definition line with no entries following \n");
           try {
            if (cur_token=="variableStep")
            {
                _processVariabledWigLine(ss);
            }
            else
            {
                _processFixedWigLine(ss);
            }
            }
            catch(uParser_missing_mandatory_values &e){
                throw e;
            }
        }
        else
            throw uParser_missing_mandatory_header()<<string_error("Missing header in wiggle file \n");

        }
}


uToken uParserWig::getNextEntry()
{
        bool foundToken=false;
        try
        {
            while(foundToken==false)
            {
                char line[4096];
                if (m_pIostream->getline(line, 4096))
                {
                    std::stringstream ss;
                    ss << line;

                    std::string cur_token;
                    ss>>cur_token;

                    if((cur_token=="variableStep")||(cur_token=="fixedStep"))
                    {
                        if (m_pIostream->eof())
                            throw end_of_file_throw()<<string_error("Badly formed wig file, track definition line with no entries following \n");
                        if (cur_token=="variableStep")
                        {
                            _processVariabledWigLine(ss);
                        }
                        else
                        {
                            _processFixedWigLine(ss);
                        }
                        //Process next line
                    }
                    else
                    {
                        //If not eof
                        int end_pos;
                        int start_pos =0;
                        float score=0.0f;
                      //  std::cerr << "Getting info " <<std::endl;
                        switch (m_Info.getStepType())
                        {
                        case wigInformation::stepType::NA:
                            throw uParser_missing_mandatory_values()<<string_error("No declaraction line in Wig, error parsing \n");
                            break;
                        case  wigInformation::stepType::FIXED:
                            if (!ss.eof())
                                throw uParser_invalid_line()<<string_error("Invalid line in file \n");

                            start_pos=m_Info.getCurPos();
                            score=stof(cur_token);
                            end_pos= start_pos+m_Info.getSpan();
                            m_Info.setCurPos(start_pos+m_Info.getStep());
                            break;
                        case   wigInformation::stepType::VARIABLE:
                            start_pos=stoi(cur_token);
                            ss >> cur_token;
                            score=stof(cur_token);
                            end_pos= start_pos+m_Info.getSpan();
                            if (!ss.eof())
                                throw uParser_invalid_line()<<string_error("Invalid line in file \n");

                            break;
                        }
                        std::stringstream token_infos;
                        token_infos << "CHR\t" << m_Info.getChrom() << "\n";
                        token_infos << "START_POS\t" << start_pos << "\n";
                        token_infos << "END_POS\t" << end_pos << "\n";
                        token_infos << "SCORE\t" << score << "\n";
                        foundToken=true;
                        return uToken(token_infos);
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
            }
        }
        catch(invalid_uToken_throw& e)
        {
            throw e;
        }
        catch(ugene_exception_base& e)
        {
           // std::cerr << "Throwing in getnextEntry uParser" <<std::endl;
            //std::cerr << fetchStringError(e) <<std::endl;
            throw e;
        }
        catch(std::exception & e)
        {
          //  std::cerr << "Throwing in getnextEntry, exception" <<std::endl;
            throw e;
        }
        std::cerr <<"Fatal error in _getNextEntryCustom(), should not reach here" <<std::endl;
        abort();
}

/** \brief Private function to process a Fixed track definition line in a wig file. Sets m_Info details
 * \return
 */
void uParserWig::_processFixedWigLine(std::stringstream & curSStream)
{
    const std::string SPANSYMBOL="span=";
    const std::string STEPSYMBOL="step=";
    const std::string CHROMSYMBOL="chrom=";
    const std::string STARTSYMBOL="start=";

    std::string chrom, strstart, strsspan;
    int curStart;
    /**< Chrom */
    curSStream >> chrom;
    if (chrom.find(CHROMSYMBOL)==std::string::npos)
        throw uParser_missing_mandatory_values()<<string_error("Missing chrom value in wig track definition");
    chrom=chrom.substr(chrom.find(CHROMSYMBOL)+CHROMSYMBOL.size());
    /**< Start */
    curSStream>>strstart;
    if (strstart.find(STARTSYMBOL)==std::string::npos)
        throw uParser_missing_mandatory_values()<<string_error("Missing start value in wig track definition");;

    strstart=strstart.substr(strstart.find(STARTSYMBOL)+STARTSYMBOL.size());
    curStart=stoi(strstart);
    /**< Step */
    std::string step;
    curSStream >>step;
    //Format definition is not sure if step is mandatory, so we will pretend it is not...
    int curStep=1;
    if (step.size())
    {
        if (step.find(STEPSYMBOL)==std::string::npos)
            throw  uParser_missing_mandatory_values()<<string_error("Missing step value in wig track definition");;;

        step=step.substr(step.find(STEPSYMBOL)+STEPSYMBOL.size());

        curStep=stoi(step);
    }
    /**< Optional Span parameter */
    int curSpan=1;
    std::string span;
    curSStream>>span;
    if(span.size())
    {
        /**< If not, fail again*/
        if (span.find(SPANSYMBOL)==std::string::npos)
            throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");

        span=span.substr(span.find(SPANSYMBOL)+SPANSYMBOL.size());
        curSpan=stoi(span);

    }
    //Nothing threw, modify values
    //Specification say's you cannot change spans
    if ((m_Info.getSpan()!=-1)&&(m_Info.getSpan()!=curSpan))
        throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file. Specification forbids changin span in dataset \n");
    m_Info.setStepType(wigInformation::stepType::FIXED);
    m_Info.setChrom(chrom);
    m_Info.setStep(curStep);
    m_Info.setSpan(curSpan);
    m_Info.setCurPos(curStart);
}

void uParserWig::_processVariabledWigLine(std::stringstream & curSStream)
{
    const std::string SPANSYMBOL="span=";
    const std::string CHROMSYMBOL="chrom=";

    std::string chrom;
    std::string span;
    int curSpan=0;
    /**< Chrom tag */
    curSStream>> chrom;
    /**< If invalid chrom header */
    if (!(chrom.size())||((chrom.find(CHROMSYMBOL)==std::string::npos)))
        throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");
    chrom=chrom.substr(chrom.find(CHROMSYMBOL)+CHROMSYMBOL.size());
    /**< Optional Span parameter */
    curSpan=1;
    curSStream >> span;
    if(span.size())
    {
        /**< If good, yay, if not, fail again*/
        if (span.find(SPANSYMBOL)==std::string::npos)
            throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");;
        //   int testpos= span.find(SPANSYMBOL);
        span=span.substr(span.find(SPANSYMBOL)+SPANSYMBOL.size());
        curSpan=stoi(span);
    }
    if ((m_Info.getSpan()!=-1)&&(m_Info.getSpan()!=curSpan))
        throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file. Specification forbids changin span in dataset \n");
    m_Info.setStepType(wigInformation::stepType::VARIABLE);
    m_Info.setChrom(chrom);
    m_Info.setSpan(curSpan);
}

DerivedParserRegister<uParserWig> uParserWig::reg("WIG");

}
