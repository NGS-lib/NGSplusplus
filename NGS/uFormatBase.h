#ifndef UFORMATBASE_H_INCLUDED
#define UFORMATBASE_H_INCLUDED

#include <string>
#include <iostream>
#include <vector>
#include "utility/utility.h"
#include "uGeneException.h"
#include "IO/uToken.h"
#include "IO/Writer/uWriter.h"

namespace NGS {
//Our basic Site for NGS format
//Very weak class as there are many differences in the functionality of derived classes.
enum class StrandDir
{
    FORWARD,REVERSE
};
enum class SplitType
{
    STRICT, IGNORE, EXTEND, ADD
};

template<class _SELF_>
class uGenericNGS
{

#define FORWARCHARD '+'
#define REVERSECHAR '-'

protected:
    std::string chr="";
    long int startPos=0;
    long int endPos=0;
    StrandDir strand=StrandDir::FORWARD;
    std::vector<float> score={};
public:
    /**< Constructor taking chromosome name, start and end */
    uGenericNGS(std::string pchr, int pstart, int pend, StrandDir pstrand=StrandDir::FORWARD ):chr(pchr),strand(pstrand)
    {
        try
        {
            setEnd(pend);
            setStart(pstart);
        }
        catch( ugene_exception_base &e)
        {
#ifdef DEBUG
            std::cerr << "Error in uGenericNGS(std::string pchr, int pstart, int pend). data is"<< pchr<< " "<< pend<<" "<< pstart << " "<<std::endl;
#endif
//            e<<generic_error(*this);
            throw e;
        }
    };

    uGenericNGS(std::string pchr, int pstart, int pend, StrandDir pstrand, float pScore ):chr(pchr),strand(pstrand)
    {
        try
        {
            setScore(pScore);
            setEnd(pend);
            setStart(pstart);
        }
        catch( ugene_exception_base &e)
        {
#ifdef DEBUG
            std::cerr << "Error in uGenericNGS(std::string pchr, int pstart, int pend, float Score). data is"<< pchr<< " "<< pend<<" "<< pstart << " "<<std::endl;
#endif
          //  e<<generic_error(*this);
            throw e;
        }
    };

    uGenericNGS(std::string pchr, int pstart, int pend,float pScore ):chr(pchr),strand(StrandDir::FORWARD)
    {
        try
        {
            setScore(pScore);
            setEnd(pend);
            setStart(pstart);
        }
        catch( ugene_exception_base &e)
        {
#ifdef DEBUG
            std::cerr << "Error in uGenericNGS(std::string pchr, int pstart, float Score). data is"<< pchr<< " "<< pend<<" "<< pstart << " "<<std::endl;
#endif
//            e<<generic_error(*this);
            throw e;
        }
    };

    uGenericNGS()
    {};
    uGenericNGS(const uToken & pToken):chr(pToken.getParam(token_param::CHR))
    {
        try
        {
            setEnd(utility::stoi(pToken.getParam(token_param::END_POS)));
            setStart( utility::stoi(pToken.getParam(token_param::START_POS)));
            /**< Default forward */
            if (pToken.isParamSet(token_param::STRAND))
                setStrand(pToken.getParam(token_param::STRAND).at(0));
            if ((pToken.isParamSet(token_param::SCORE))&&(pToken.getParam(token_param::SCORE)!="." ) )  {
                    setScore(utility::stof (pToken.getParam(token_param::SCORE) ) );
                }
        }
        catch(ugene_exception_base &e)
        {
#ifdef DEBUG
            std::cerr << "Error in uGenericNGS(uToken)." <<std::endl;
#endif
            e<<string_error("Error in uGenericNGS(uToken)." );
            throw e;
        }

    };

    virtual ~uGenericNGS() {};
    /**< End Constructor/Destructor */

    /**< Write to Bed functions, to be deprecated*/
    virtual void writeBedToOuput(std::ostream &out, bool endLine) const;
    virtual void writeBedToOuput(std::ostream &out) const
    {
        writeBedToOuput(out, true);
    };
     virtual void writeToOutput(uWriter& pWriter) const;
     virtual uToken getToken()const;

    /**< Get /Set */
    std::string getChr() const
    {
        return chr;
    };
    void setChr(std::string ourchr)
    {
        if (ourchr.size()==0)
            throw param_throw() << string_error("Throwing in setChr, ID must be of size > 0");
        chr=ourchr;
    };

    StrandDir getStrand() const
    {
        return strand;
    };

    void setStrand(char pStrand)
    {
        try
        {
            if (pStrand==REVERSECHAR)
                strand=StrandDir::REVERSE;
            else if (pStrand==FORWARCHARD)
                strand=StrandDir::FORWARD;
            else
                throw param_throw()<<string_error("Failling in setStrand(char), invalid character");
        }
        catch(param_throw &e)
        {
            //elem_throw e;
            throw e;
        }
    };
    void setStrand(StrandDir pStrand)
    {
        strand=pStrand;
    };

    bool isReverse() const
    {
        if (strand==StrandDir::REVERSE)
            return true;
        else
            return false;
    }

    void setStart(int ourStart)
    {
        try
        {
            if (!((ourStart<=getEnd())&&(ourStart>=0)))
                throw param_throw()<<string_error("Failed in setStart, ourStart is smalled then end or under 0, start is "+utility::to_string(ourStart)+ " end is "+ utility::to_string(getEnd()) +"\n");
            startPos=ourStart;

        }
        catch(param_throw &e)
        {
#ifdef DEBUG
            std::cerr << "throwing in setStart" <<std::endl;
#endif
            throw e;
        }
    };
    void setEnd(int ourEnd)
    {
        try
        {
            if (!((ourEnd>=getStart())&&(ourEnd>=0)))
                throw param_throw()<<string_error("throwing in setEnd(), start at "+utility::to_string((int)getStart())+ " end is "+ utility::to_string(ourEnd) +"\n");
            endPos=ourEnd;
        }
        catch(param_throw & e)
        {
#ifdef DEBUG
            std::cerr << "throwing in setEnd" <<std::endl;
#endif
            throw e;
        }
    };

    void setStartEnd(long int ourStart, long int ourEnd)
    {
        try
        {
        if ((ourStart>0) && (ourStart<=ourEnd))
            {
            startPos=ourStart;
            endPos=ourEnd;
            }
        else
            {
            throw param_throw() <<string_error("Throwing in set StartEnd,  ourStart="+utility::to_string(ourStart)+" ourEnd="+utility::to_string(ourEnd)+"\n" );
            }
        }
        catch(param_throw &e)
        {
#ifdef DEBUG
            std::cerr << "throwing in setStartEnd" <<std::endl;
#endif
            throw e;
        }
    };

    long int getStart() const
    {
        return startPos;
    };
    long int getEnd() const
    {
        return endPos;
    };
    long int getLenght() const
    {
        /**< 0 based coordinates, so N - N  is a legal fragment covering a single nucleotide at position N */
        return (endPos-startPos+1);
    };


    /**< Strictly for debugging */
    virtual void debugElem() const
    {
        using namespace utility;
        stringTocerr("Outputting elemn data");
        stringTocerr("Chrom "+getChr());
        stringTocerr("Start "+utility::to_string((int)getStart()));
        stringTocerr("End " +utility::to_string((int)getEnd()));
    }
    /**<  Divide our region into a certain number of subregions */

    std::vector<_SELF_> divideIntoBinofSize(const int N, const SplitType type=SplitType::STRICT);
    std::vector<_SELF_> divideIntoNBin(const int N, const SplitType ptype=SplitType::STRICT);

    //If we want to Change the dimensions of our site
    void extendSite(int extend);
    void extendSite(int extendLeft, int extendRight);

    void trimSites(int trim);
    void trimSites(int trimLeft,int trimRight);

    /**< Should this be there? */
    bool doesOverlap(_SELF_ other,OverlapType type=OverlapType::OVERLAP_PARTIAL) const;

    float getScore(int p_Pos) const;
    float getScore()const{return getScore(0);};
    int getScoreCount() const { return score.size();};

    void setScore(float p_score, int p_Pos);
    void setScore(float ourscore) {setScore(ourscore,0);}

};

/** \brief Increase size of the element. Coordinates can go no lower then 0,
 *
 * \param extendLeft int : Left shift size, must be +
 * \param extendRight int : Right shift size, must be +
 * \return void
 *
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::extendSite(int extendLeft, int extendRight)
{
    try
    {
        if((extendLeft<0)||(extendRight<0))
            throw param_throw()<< string_error("INIT throwing from extendSite("+utility::to_string(extendLeft)+","+utility::to_string(extendRight)+"), param < 0 \n"  );;

        int start=(getStart()-extendLeft);
        if (start < 0)
            start=0;

        setStart(start);
        setEnd(getEnd()+extendRight);
    }
    catch(param_throw & e )
    {
        elem_throw e;
        std::string * trace;

     if ( (trace=(boost::get_error_info<string_error>(e))) )
        e << string_error(*trace+"Catching and re-throwing from extendSite("+utility::to_string(extendLeft)+","+utility::to_string(extendRight)+")\n");
     else
         e << string_error("Catching and re-throwing from extendSite("+utility::to_string(extendLeft)+","+utility::to_string(extendRight)+")\n");
        throw(e);
        return;
    }
}
/** \brief Idem as above, but equal shift.
 *
 * \param extend int : Size of shift on both sides
 * \return void
 *
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::extendSite(int extend)
{
    try
    {
        this->extendSite(extend,extend);
    }
    catch(param_throw & e)
    {
        throw e;
    }

}

/** \brief Diminish our element by reducing from either or both sides.
 *
 * \param trimLeft int : Size to trim from left coord
 * \param trimRight int : Size to trim from right cood
 * \pre Cannot trim a negative. Cannot trim so that start and end swap
 * \
 * \return void
 *
 */
 template <class _SELF_>
void uGenericNGS<_SELF_>::trimSites(int trimLeft, int trimRight)
{
    /**< Validate input */
    try
    {
        if ((trimLeft<0)||(trimRight<0)||(trimLeft+trimRight>this->getLenght()))
            throw param_throw()<< string_error("PARAMERROR, throwing from trimSites("+utility::to_string(trimLeft)+","+utility::to_string(trimRight)+"), param < 0 \n"  );

        this->startPos=(this->startPos+trimLeft);
        this->endPos=(this->endPos-trimRight);
    }
    catch (param_throw & err)
    {
        #ifdef DEBUG
           std::cerr << "Throwing in trimSites(int, int)"<<std::endl;
        #endif
        throw(err);
    }


}

/** \brief Diminish our element by an equal amount from both ends
 *
 * \param trim int : Size to trim from each side
 * \return void
 *
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::trimSites(int trim)
{
    try {
        this->trimSites(trim,trim);
    }
    catch(param_throw & e)
    {
        throw e;
    }

}

/** \brief Verify if two sites overlap by at least 1 bp
 *
 * \param other uGenericNGS : Other site to validate
 * \return bool : True if overlap.
 *
 */
template <class _SELF_>
bool uGenericNGS<_SELF_>::doesOverlap(_SELF_ other, OverlapType type) const
{
    bool returnb=false;
    if (getChr()==other.getChr())
        returnb=utility::isOverlap(this->startPos, this->endPos, other.startPos,other.endPos,type);

    return returnb;
}

/** \brief Deprecated code, output in Bed format, Endline can be optionally be delegated
 *  \brief Will be removed next update
 *
 * \param out std::ostream& : Our output steram
 * \param endLine bool : If true, write endline
 * \return void
 *
 */
 template <class _SELF_>
void uGenericNGS<_SELF_>::writeBedToOuput(std::ostream &out, bool endLine) const
{

    out << chr << "\t" << startPos << "\t" << endPos;

    if (endLine)
        out <<std::endl;
}

/** \brief Write contig to writer.
 *
 * \param pWriter uWriter& Previously declared writer, passe by param
 * \return void
 *
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::writeToOutput(uWriter& pWriter) const
{
    pWriter.writeToken(this->getToken());
}



/** \brief Divide the contig into N contigs of equal size
 *
 * \param N int : Number of entries to create from the entry
 * \return vector<uGenericNGS> List of created contigs
 *
 */
  template <class _SELF_>
std::vector<_SELF_> uGenericNGS<_SELF_>::divideIntoNBin(int N,SplitType ptype)
{
    std::vector<_SELF_> returnVec;
    int leftover = getLenght()%N;
    int binSize= getLenght()/N;
    try
    {
        /**< If NB bins is greater then BP */
        if (this->getLenght()<N)
           throw param_throw() << string_error("Asking for more bins then lenght of Elem /n");

        /**< If Strict and we cannot exactly fit our bins, fail */
        if ((ptype==SplitType::STRICT)&&(leftover!=0))
            throw param_throw()<< string_error( "STRICT parameter not respected /n");

        /**< Once we pass the check, generate our basic bin */
        int curStart=getStart();
        for (int i=0; i<(N); i++)
        {
            returnVec.push_back(  _SELF_(getChr(),curStart, (curStart+(binSize-1))  ));
            curStart+=binSize;
        }
        /**< Then we deal with any leftovers according to policy */
        /**< If IGNORE, we do nothing */
        if(leftover!=0)
        {
            switch(ptype)
            {
                /**< Extend the final region to fill the gap */
            case SplitType::EXTEND:
            {
                returnVec.at(N-1).setEnd(curStart+(leftover-1));
                break;
            }
            /**< Add an extra region for the leftover */
            case SplitType::ADD:
            {
                returnVec.push_back(  _SELF_(getChr(),curStart, (curStart+(leftover-1))  ));
                break;
            }
            case SplitType::IGNORE:
            {
            break;
            }
            case SplitType::STRICT:
            {
                throw param_throw()<<string_error("Invalid trace in divideIntoNBin, STRICT should have been validated earlier");
            break;
            }

            }
        }
    }
    catch(param_throw &e)
    {
       std::string ourtrace="Failed in divideIntoNBin()";
        if(  std::string const * trace=boost::get_error_info<string_error>(e) )
            ourtrace+=*trace;

        #ifdef DEBUG
                std::cerr <<"Failed in divideIntoNBin()"<<std::endl;
        #endif
        e << string_error(ourtrace);
        throw e;
    }
    catch(std::exception &e ){
        #ifdef DEBUG
                std::cerr <<"Failed in divideIntoNBin()"<<std::endl;
        #endif
        throw e ;
    }

    return returnVec;
}

/** \brief Divide the contig into contigs of size N
 *
 * \param N int : Size of bins to create.
 * \param const SplitType: Determines how you managed leftover
 * \return vector<uGenericNGS> List of returned contigs
 *
 */
 template <class _SELF_>
std::vector<_SELF_> uGenericNGS<_SELF_>::divideIntoBinofSize(const int N, const SplitType type)
{
    std::vector<_SELF_> returnVec;
    try
    {
        int binSize= N;
        int nbBin= getLenght()/N;
        int leftover = getLenght()%binSize;

        if (getLenght()<N)
            throw param_throw() << string_error("Asking for more bins then lenght of Elem /n");

        /**< If Strict and we cannot exactly fit our bins, fail */
        if ((type==SplitType::STRICT)&&(leftover!=0))
            throw param_throw()<< string_error( "STRICT parameter not respected /n");

        /**< Once we pass the check, generate our basic bin */

        int curStart=getStart();
        for (int i=0; i<nbBin; i++)
        {
            returnVec.push_back(  uGenericNGS(getChr(),curStart, (curStart+(binSize-1))  ));
            curStart+=binSize;
        }

        /**< Then we deal with any leftovers according to policy */
        /**< If IGNORE, we do nothing */
        if(leftover!=0)
            switch(type)
            {
                /**< Extend the final region to fill the gap */
            case SplitType::EXTEND:
            {
                returnVec.at(nbBin-1).setEnd(curStart+(leftover-1));
                break;
            }
            /**< Add an extra region for the leftover */
            case SplitType::ADD:
            {
                returnVec.push_back(  uGenericNGS(getChr(),curStart, (curStart+(leftover-1))  ));
                break;
            }
             case SplitType::IGNORE:
            {
            break;
            }
            case SplitType::STRICT:
            {
                throw param_throw()<<string_error("Invalid trace in divideIntoNBin, STRICT should have been validated earlier");
            break;
            }
            }
    }
    catch(param_throw &e)
    {
         std::string ourtrace="Failed in divideIntoNBin() /n";
        if(  std::string const * trace=boost::get_error_info<string_error>(e) )
            ourtrace+=*trace;

        #ifdef DEBUG
                std::cerr <<"Failed in divideIntoBinofSize()"<<std::endl;
        #endif
        e << string_error(ourtrace);
        throw e;
    }
    catch(std::exception & e)
    {
         #ifdef DEBUG
                  std::cerr <<"Failed in divideIntoBinofSize()"<<std::endl;
        #endif
        throw e;
    }

    return returnVec;
}
/**< A score is an arbitray value set that can be used later */
/** \brief Set the score of a contig. Note that this involves resizing the vector, so settting arbitrarily large score counts can bust your memory
 * \param float p_score: the score to set
 * \param int p_Pos: the count of the score to set
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::setScore(float p_score, int p_Pos)
{
    try {
    if (p_Pos>= ((int)score.size()))
        score.resize(p_Pos+1);
    score.at(p_Pos)=p_score;
    }
    catch(std::exception &e){throw e;}
}

/** \brief Fetch the score at a specific region, return infinity if not set
 * \param int p_Pos: the position where the score should be fetched
 */
template <class _SELF_>
float uGenericNGS<_SELF_>::getScore(int p_Pos) const
{
    if (p_Pos>=((int)score.size()))
       return std::numeric_limits<float>::infinity();
    return score.at(p_Pos);
}


 /** \brief Get the parser Token associated with our Read/Write class
  *
  * \return uToken The token returned
  *
  */
template <class _SELF_>
 uToken uGenericNGS<_SELF_>::getToken() const
 {
    std::stringstream ss;
    ss << "CHR\t"<<this->getChr()<<"\nSTART_POS\t"<<this->getStart()<<"\n" << "END_POS\t"<<this->getEnd()<<"\n";
    if (getScore()!=std::numeric_limits<float>::infinity())
       ss << "SCORE\t"<<this->getScore()<<"\n";

    return uToken(ss);
 }



} // End of namespace NGS

/**< Legacy code, will return a uGenericNGS from  a tab delimited code. Deprecated, will be removed */
namespace factory
{
    template<class _SELF_>
    inline static _SELF_ makeNGSfromTabString(std::string tabString)
    {

        utility::Tokenizer tabLine(tabString);

        std::string chrm;
        int start=0, end=0;
        try
        {
            tabLine.NextToken();
            chrm = tabLine.GetToken();
            tabLine.NextToken();
            start = utility::stoi(tabLine.GetToken());
            tabLine.NextToken();
            end = utility::stoi(tabLine.GetToken());

            return _SELF_(chrm,start,end);
        }
        catch(...)
        {
            throw;
        }

    }
}

#endif // UFORMATBASE_H_INCLUDED
