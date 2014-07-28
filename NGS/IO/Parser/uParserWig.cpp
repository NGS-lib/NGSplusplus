// ***************************************************************************
// uParserWig.cpp (c) 2013
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


#include "uParserWig.h"
#include "../../utility/utility.h"
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
        std::string line;
    if (std::getline(*m_pIostream,line))
    {
		utility::GetTokens(m_tokens,line);
        std::string cur_token=m_tokens.at(0);

		/**< Track line, get next one, then fail if not valid */
		if(cur_token=="track"){
			if (std::getline(*m_pIostream,line))
				{
					utility::GetTokens(m_tokens,line);
					cur_token=m_tokens.at(0);
				}
			else
				throw uParser_missing_mandatory_header()<<string_error("Missing header in wiggle file \n");
		}
        if((cur_token=="variableStep")||(cur_token=="fixedStep"))
        {
            if (m_pIostream->eof())
                throw end_of_file_throw()<<string_error("Badly formed wig file, track definition line with no entries following \n");
           try {
                if (cur_token=="variableStep")
                {
                    _processVariabledWigLine(m_tokens);
                }
                else
                {
                    _processFixedWigLine(m_tokens);
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
void uParserWig::init(std::istream* stream, bool header )
{
    m_pIostream = stream;
    m_dynamicStream = false;
     std::string line;
    if (std::getline(*m_pIostream,line))
    {
     //  std::stringstream ss;
//        ss << line;

     //   ss>>cur_token;
		utility::GetTokens(m_tokens,line);
        std::string cur_token=m_tokens.at(0);

		if(cur_token=="track"){
			if (std::getline(*m_pIostream,line))
				{
					utility::GetTokens(m_tokens,line);
					cur_token=m_tokens.at(0);
				}
			else
				throw uParser_missing_mandatory_header()<<string_error("Missing header in wiggle file \n");
		}

        if((cur_token=="variableStep")||(cur_token=="fixedStep"))
        {
            if (m_pIostream->eof())
                throw end_of_file_throw()<<string_error("Badly formed wig file, track definition line with no entries following \n");
           try {
            if (cur_token=="variableStep")
            {
                _processVariabledWigLine(m_tokens);
            }
            else
            {
                _processFixedWigLine(m_tokens);
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
                std::string line;
                if (std::getline(*m_pIostream,line))
                {
                     m_rawString=line;
                    utility::GetTokens(m_tokens,line);
                    if((m_tokens.at(0)=="variableStep")||(m_tokens.at(0)=="fixedStep"))
                    {
                        if (m_pIostream->eof())
                            throw end_of_file_throw()<<string_error("Badly formed wig file, track definition line with no entries following \n");
                        if (m_tokens.at(0)=="variableStep")
                        {
                            _processVariabledWigLine(m_tokens);
                        }
                        else
                        {
                            _processFixedWigLine(m_tokens);
                        }
                        //Process next line
                    }
                    else
                    {
                        //If not eof
                        int end_pos=0;
                        int start_pos =0;
                        float score;
                        int scorePos=0;
                        switch (m_Info.getStepType())
                        {
                        case wigInformation::stepType::NA:
                            throw uParser_missing_mandatory_values()<<string_error("No declaraction line in Wig, error parsing \n");
                            break;
                        case  wigInformation::stepType::FIXED:
                            //if (!ss.eof())
                            if(m_tokens.size()!=1)
                                throw uParser_invalid_line()<<string_error("Invalid line in file \n");

                            start_pos=m_Info.getCurPos();
                            scorePos=0;
                            score=utility::stof(m_tokens.at(scorePos));
                            end_pos= start_pos+m_Info.getSpan();
                            m_Info.setCurPos(start_pos+m_Info.getStep());
                            break;
                        case   wigInformation::stepType::VARIABLE:
                            start_pos=utility::stoi(m_tokens.at(0));
							scorePos=1;
                            score=utility::stof(m_tokens.at(scorePos));
                            end_pos= start_pos+m_Info.getSpan();

                            if(m_tokens.size()>2)
                                throw uParser_invalid_line()<<string_error("Invalid line in file \n");

                            break;
                        }
                        uToken ourToken;
                        ourToken._setParamNoValidate(token_param::CHR,m_Info.getChrom());
						ourToken._setParamNoValidate(token_param::START_POS,utility::to_string(start_pos));
						ourToken._setParamNoValidate(token_param::END_POS,utility::to_string(end_pos));
						ourToken._setParamNoValidate(token_param::SCORE,utility::to_string(score));
                        foundToken=true;

                        return ourToken;
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
            throw e;
        }
        catch(std::exception & e)
        {
          throw uParser_exception_base()<<string_error(e.what());
        }
        std::cerr <<"Fatal error in _getNextEntryCustom(), should not reach here" <<std::endl;
        abort();
}

/** \brief Private function to process a Fixed track definition line in a wig file. Sets m_Info details
 * \return
 */
void uParserWig::_processFixedWigLine(std::vector<std::string> & curSStream)
{
    const std::string SPANSYMBOL="span=";
    const std::string STEPSYMBOL="step=";
    const std::string CHROMSYMBOL="chrom=";
    const std::string STARTSYMBOL="start=";

    std::string chrom, strstart, strsspan;
    int curStart;
    /**< Chrom */
     chrom=curSStream.at(1);
    if (chrom.find(CHROMSYMBOL)==std::string::npos)
        throw uParser_missing_mandatory_values()<<string_error("Missing chrom value in wig track definition");
    chrom=chrom.substr(chrom.find(CHROMSYMBOL)+CHROMSYMBOL.size());
    /**< Start */
    strstart=curSStream.at(2);
    if (strstart.find(STARTSYMBOL)==std::string::npos)
        throw uParser_missing_mandatory_values()<<string_error("Missing start value in wig track definition");;

    strstart=strstart.substr(strstart.find(STARTSYMBOL)+STARTSYMBOL.size());
    curStart=utility::stoi(strstart);
    /**< Step */
    std::string step;
    if ( curSStream.size()>3)
    	step = curSStream.at(3);
    //Format definition is not sure if step is mandatory, so we will pretend it is not...
    int curStep=1;
    if (step.size())
    {
        if (step.find(STEPSYMBOL)==std::string::npos)
            throw  uParser_missing_mandatory_values()<<string_error("Missing step value in wig track definition");

        step=step.substr(step.find(STEPSYMBOL)+STEPSYMBOL.size());
        curStep=utility::stoi(step);
    }
    /**< Optional Span parameter */
    int curSpan=1;
    std::string span;
    if ( curSStream.size()>4)
    	span = curSStream.at(4);
    if(span.size())
    {
        /**< If not, fail again*/
        if (span.find(SPANSYMBOL)==std::string::npos)
            throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");

        span=span.substr(span.find(SPANSYMBOL)+SPANSYMBOL.size());
        curSpan=utility::stoi(span);

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

void uParserWig::_processVariabledWigLine(std::vector<std::string> & curSStream)
{
    const std::string SPANSYMBOL="span=";
    const std::string CHROMSYMBOL="chrom=";

    std::string chrom;
    std::string span;
    int curSpan=0;
    /**< Chrom tag */
	chrom=curSStream.at(1);
    /**< If invalid chrom header */
    if (!(chrom.size())||((chrom.find(CHROMSYMBOL)==std::string::npos)))
        throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");
    chrom=chrom.substr(chrom.find(CHROMSYMBOL)+CHROMSYMBOL.size());
    /**< Optional Span parameter */
    curSpan=1;
    if (curSStream.size()>1){
        span=curSStream.at(2);}
    if(span.size())
    {
        /**< If good, yay, if not, fail again*/
        if (span.find(SPANSYMBOL)==std::string::npos)
            throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file, failling \n");;
        //   int testpos= span.find(SPANSYMBOL);
        span=span.substr(span.find(SPANSYMBOL)+SPANSYMBOL.size());
        curSpan=utility::stoi(span);
    }
    if ((m_Info.getSpan()!=-1)&&(m_Info.getSpan()!=curSpan))
        throw uParser_missing_mandatory_values()<<string_error("invalid Track definition line in wig file. Specification forbids changin span in dataset \n");
    m_Info.setStepType(wigInformation::stepType::VARIABLE);
    m_Info.setChrom(chrom);
    m_Info.setSpan(curSpan);
}


//DerivedParserRegister<uParserWig> uParserWig::reg("WIG");

}
