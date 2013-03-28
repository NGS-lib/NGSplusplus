#ifndef UFORMATBASE_H_INCLUDED
#define UFORMATBASE_H_INCLUDED

#include <string>
#include <iostream>
#include <vector>
#include "utility/utility.h"
#include "uGeneException.h"
#include "IO/uToken.h"
#include "IO/Writer/uWriter.h"

namespace NGS
{
/**< Our basic Site for NGS format */
/**< Very weak class as there are many differences in the functionality of derived classes. */

enum class StrandDir: bool /*!< Possible strand directions */
{
    FORWARD,REVERSE
};

enum class SplitType /*!< Used with division functions */
{
    STRICT, IGNORE, EXTEND, ADD
};

template<class _SELF_>
class uGenericNGS
{

#define FORWARDCHAR '+'
#define REVERSECHAR '-'

protected:
    std::string m_chr=""; /*!< Name of the scaffold/chromosome  */
    long int m_startPos=0; /*!< Begining genomic position of the element  */
    long int m_endPos=0; /*!< End genomic position of the element  */
    StrandDir m_strand=StrandDir::FORWARD; /*!< Strand orientation of the element  */
    std::vector<float> m_score= {}; /*!< Values potentially associated with the element */

public:

    uGenericNGS(std::string pChr, long int pStart, long int pEnd, StrandDir pStrand=StrandDir::FORWARD );
    uGenericNGS(std::string pChr, long int pStart, long int pEnd, StrandDir pStrand, float pScore );
    uGenericNGS(std::string pChr, long int pStart, long int pEnd,float pScore );
    uGenericNGS(){};
    uGenericNGS(const uToken & pToken);

    virtual ~uGenericNGS() {};

    /**< End Constructor/Destructor */

    virtual void writeToOutput(uWriter& pWriter) const;
    virtual uToken createToken()const;
    virtual void print(std::ostream &pOut)const;

    /**< Get/Set */

    std::string getChr() const;
    void setChr(std::string ourm_chr);

    StrandDir getStrand() const;
    void setStrand(char pStrand);
    void setStrand(StrandDir pStrand);

    bool isReverse() const;

    void setStart(long int ourStart);
    void setEnd(long int ourEnd);
    void setStartEnd(long int ourStart, long int ourEnd);

    long int getStart() const;
    long int getEnd() const;


    long int getLenght() const;

    //  virtual bool isEqual(const  _SELF_ & pCompared)const =0;
    /**<  Divide our region into a certain number of subregions */
    std::vector<_SELF_> divideIntoBinofSize(const int N, const SplitType type=SplitType::STRICT);
    std::vector<_SELF_> divideIntoNBin(const int N, const SplitType ptype=SplitType::STRICT);

    //If we want to Change the dimensions of our site
    void extendSite(long int extend);
    void extendSite(long int extendLeft, long int extendRight);
    void trimSite(long int trim);
    void trimSite(long int trimLeft,long int trimRight);

    /**< Should this be there? */
    bool doesOverlap(_SELF_ other,OverlapType type=OverlapType::OVERLAP_PARTIAL) const;

    float getScore(int p_Pos) const;
    float getScore()const
    {
        return getScore(0);
    };
    int getScoreCount() const
    {
        return m_score.size();
    };

    void setScore(float p_score, int p_Pos);
    void setScore(float ourscore)
    {
        setScore(ourscore,0);
    }
    std::vector<float> getScoreVector()const
    {
        return m_score;
    };
    void setScoreVector(std::vector<float> p_Score)
    {
        m_score=std::move(p_Score);
    };

    _SELF_ returnOverlapping(const _SELF_ &)const;
    _SELF_ returnMerge(const _SELF_ &)const;

};


/** \brief If overlapping, return the parts overlapping, if not raise exception. Suggested to always test overlap before calling
 *
 * \param _SELF_ element to compare overlap with
 * \exception param_throw() : Throw if the param element does not overlap.
 * \return _SELF_ The element returned
 *
 */
template <class _SELF_>
_SELF_ uGenericNGS<_SELF_>::returnOverlapping(const _SELF_ & otherItem)const
{
    if (utility::isOverlap(this->getStart(),this->getEnd(), otherItem.getStart(), otherItem.getEnd()))
    {
        long long curStart =this->getStart(), curEnd=this->getEnd();
        StrandDir dir=StrandDir::FORWARD;
        /**< Only keep strand direction if both items agree, otherwise set to FORWARD */
        if (this->getStrand()==otherItem.getStrand())
            dir= this->getStrand();
        if (this->getStart()<otherItem.getStart())
            curStart=otherItem.getStart();
        if (this->getEnd()>otherItem.getEnd())
            curEnd=otherItem.getEnd();

        return (_SELF_(this->getChr(),curStart,curEnd,dir));
    }
    else
        throw param_throw();
}

/** \brief Return an element with the boundaries merged between this and comparison.
 *
 * \param otherItem const _SELF_& : Element to merge with.
 * \exception param_throw : When element do not overlap.
 * \return _SELF_ new element with boundaries merge.
 *
 */
 template <class _SELF_>
_SELF_ uGenericNGS<_SELF_>::returnMerge(const _SELF_ & otherItem)const
{
    if (utility::isOverlap(this->getStart(),this->getEnd(), otherItem.getStart(), otherItem.getEnd()))
    {
        long long curStart =this->getStart(), curEnd=this->getEnd();
        StrandDir dir=StrandDir::FORWARD;
        /**< Only keep strand direction if both items agree, otherwise set to FORWARD */
        if (this->getStrand()==otherItem.getStrand())
            dir= this->getStrand();
        if (this->getStart()>otherItem.getStart())
            curStart=otherItem.getStart();
        if (this->getEnd()<otherItem.getEnd())
            curEnd=otherItem.getEnd();
        return (_SELF_(this->getChr(),curStart,curEnd,dir));
    }
    else
        throw param_throw();
}
/** \brief Constructor taking chromosome name, start, end and optionnaly a strand.
* \param std::string pChr: The name of the chromosome.
* \param int pStart: Starting position, must be greater than 0.
* \param int pEnd: Ending position, must be greater or equal to pStart.
* \param StrandDir pStrand: The strand in enum class StranDir format. Must be either FOWARD or REVERSE.
* \exception ugene_exception_base: When the starting position and ending position are incorrect.
*/
template <class _SELF_>
uGenericNGS<_SELF_>::uGenericNGS(std::string pChr, long int pStart, long  int pEnd, StrandDir pStrand):m_chr(pChr),m_strand(pStrand)
{
    try
    {
        setEnd(pEnd);
        setStart(pStart);
    }
    catch( ugene_exception_base &e)
    {
#ifdef DEBUG
        std::cerr << "Error in uGenericNGS(std::string pChr, int pStart, int pEnd). data is"<< pChr<< " "
                  << pEnd<<" "<< pStart << " "<<'\n';
#endif
        throw e;
    }
}

/** \brief Constructor taking chromosome name, start position, end position, strand and score.
* \param std::string pChr: The name of the chromosome.
* \param int pStart: Starting position, must be greater than 0.
* \param int pEnd: Ending position, must be greater or equal to pStart.
* \param StrandDir pStrand: The strand in enum class StranDir format. Must be either FOWARD or REVERSE.
* \param float pScore: The score associated with the current entry.
* \exception ugene_exception_base: When the starting position and ending position are incorrect.
*/
template <class _SELF_>
uGenericNGS<_SELF_>::uGenericNGS(std::string pChr, long  int pStart, long  int pEnd, StrandDir pStrand, float pScore ):m_chr(pChr),m_strand(pStrand)
{
    try
    {
        setScore(pScore);
        setEnd(pEnd);
        setStart(pStart);
    }
    catch(ugene_exception_base &e)
    {
#ifdef DEBUG
        std::cerr << "Error in uGenericNGS(std::string pChr, int pStart, int pEnd, float Score). data is"
                  << pChr<< " "<< pEnd<<" "<< pStart << " "<<std::endl;
#endif
        throw e;
    }
}

/** \brief Constructor taking chromosome name, start position, end position and score.
* \param std::string pChr: The name of the chromosome.
* \param int pStart: Starting position, must be greater than 0.
* \param int pEnd: Ending position, must be greater or equal to pStart.
* \param float pScore: The score associated with the current entry.
* \exception ugene_exception_base: When the starting position and ending position are incorrect.
*/
template <class _SELF_>
uGenericNGS<_SELF_>::uGenericNGS(std::string pChr, long  int pStart, long  int pEnd,float pScore ):m_chr(pChr),m_strand(StrandDir::FORWARD)
{
    try
    {
        setScore(pScore);
        setEnd(pEnd);
        setStart(pStart);
    }
    catch( ugene_exception_base &e)
    {
#ifdef DEBUG
        std::cerr << "Error in uGenericNGS(std::string pChr, int pStart, float Score). data is"
                  << pChr<< " "<< pEnd<<" "<< pStart << " "<<std::endl;
#endif
        throw e;
    }
}

/** \brief Constructor taking a token.
   * \param const uToken & pToken: A token containing all the infos for the current region. See uToken class for more details.
   * \exception ugene_exception_base: When the values of the token are not valid.
   */
template <class _SELF_>
uGenericNGS<_SELF_>::uGenericNGS(const uToken & pToken):m_chr(pToken.getParam(token_param::CHR))
{
    try
    {
        setEnd(utility::stoi(pToken.getParam(token_param::END_POS)));
        setStart( utility::stoi(pToken.getParam(token_param::START_POS)));
        /**< Default forward */
        if (pToken.isParamSet(token_param::STRAND))
        {
            setStrand(pToken.getParam(token_param::STRAND).at(0));
        }
        if ((pToken.isParamSet(token_param::SCORE))&&(pToken.getParam(token_param::SCORE)!="." ) )
        {
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
}


/** \brief Get the chromosome of the current entry.
 * \return A string corresponding to the name of the chr. Returns an empty string if the chr value is not set.
 */
template <class _SELF_>
std::string uGenericNGS<_SELF_>::getChr() const
{
    return m_chr;
}

/** \brief Set the name of the chromosome for the current entry.
 * \return void
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::setChr(std::string ourm_chr)
{
    if (ourm_chr.size()==0)
    {
        throw param_throw() << string_error("Throwing in setChr, ID must be of size > 0");
    }
    m_chr=ourm_chr;
}


/** \brief Set the name of the chromosome for the current entry.
 * \return A StrandDir corresponding to the strand for the current entry. The default value is FORWARD.
 */
template <class _SELF_>
StrandDir uGenericNGS<_SELF_>::getStrand() const
{
    return m_strand;
}

/** \brief Set the strand for the current entry with a char.
 * \param char pStrand: The value for the strand. Must be either '+' or '-'.
 * \exception param_throw: When the user gave an incorrect value for the strand.
 * \return void
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::setStrand(char pStrand)
{
    try
    {
        if (pStrand==REVERSECHAR)
        {
            m_strand=StrandDir::REVERSE;
        }
        else if (pStrand==FORWARDCHAR)
        {
            m_strand=StrandDir::FORWARD;
        }
        else
        {
            throw param_throw()<<string_error("Failling in setStrand(char), invalid character");
        }
    }
    catch(param_throw &e)
    {
        throw e;
    }
}

/** \brief Set the strand for the current entry with a StrandDir.
 * \param StrandDir pStrand: The value for the strand. Must be either StrandDir::FORWARD or StrandDir::REVERSE.
 * \return void
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::setStrand(StrandDir pStrand)
{
    m_strand=pStrand;
}

/** \brief Check if the orientation of the strand for the current entry is reverse.
 * \return bool: true if the strand is StrandDir::REVERSE. false if the strand is StrandDir::FORWARD.
 */
template <class _SELF_>
bool uGenericNGS<_SELF_>::isReverse() const
{
    if (m_strand==StrandDir::REVERSE)
    {
        return true;
    }

    return false;
}

/** \brief Set the start position for the current entry.
  * \param int ourStart: The start position we wish to set.
  * \exception param_throw: If the start position is greater than the end position or if the start position is smaller than zero.
  * \return void
  */
template <class _SELF_>
void uGenericNGS<_SELF_>::setStart(long int ourStart)
{
    try
    {
        if (!((ourStart<=getEnd())&&(ourStart>=0)))
        {
            throw param_throw()<<string_error("Failed in setStart, ourStart is smalled then end or under 0, start is "
                                              +utility::to_string(ourStart)+ " end is "+ utility::to_string(getEnd()) +"\n");
        }
        m_startPos=ourStart;
    }
    catch(param_throw &e)
    {
#ifdef DEBUG
        std::cerr << "throwing in setStart" <<std::endl;
#endif
        throw e;
    }
}

/** \brief Set the end position for the current entry.
 * \param int ourEnd: The end position we wish to set.
 * \exception param_throw: If the end position is smaller than the start position or if the end position is smaller than zero.
 * \return void
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::setEnd(long int ourEnd)
{
    try
    {
        if (!((ourEnd>=getStart())&&(ourEnd>=0)))
        {
            throw param_throw()<<string_error("throwing in setEnd(), start at "
                                              +utility::to_string((int)getStart())+ " end is "+ utility::to_string(ourEnd) +"\n");
        }
        m_endPos=ourEnd;
    }
    catch(param_throw & e)
    {
#ifdef DEBUG
        std::cerr << "throwing in setEnd" <<std::endl;
#endif
        throw e;
    }
}

/** \brief Set the start and end position for the current entry.
 * \param int ourStart: The start position we wish to set.
 * \param int ourEnd: The end position we wish to set.
 * \exception param_throw: If the end position is smaller than the start position or if the start and end position are smaller than zero.
 * \return void
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::setStartEnd(long int ourStart, long int ourEnd)
{
    try
    {
        if ((ourStart>0) && (ourStart<=ourEnd))
        {
            m_startPos=ourStart;
            m_endPos=ourEnd;
        }
        else
        {
            throw param_throw() <<string_error("Throwing in set StartEnd,  ourStart="
                                               +utility::to_string(ourStart)+" ourEnd="+utility::to_string(ourEnd)+"\n" );
        }
    }
    catch(param_throw &e)
    {
#ifdef DEBUG
        std::cerr << "throwing in setStartEnd" <<std::endl;
#endif
        throw e;
    }
}


/** \brief Get the start position of the current entry.
 * \return long int: the start position of the current entry. Default value is 0.
 */
template <class _SELF_>
long int uGenericNGS<_SELF_>::getStart() const
{
    return m_startPos;
}

/** \brief Get the end position of the current entry.
 * \return long int: the end position of the current entry. Default value is 0.
 */
template <class _SELF_>
long int uGenericNGS<_SELF_>::getEnd() const
{
    return m_endPos;
}

/** \brief Get the length of the current entry.
 * \return long int: the difference between the ending position and the starting position.
 */
template <class _SELF_>
long int uGenericNGS<_SELF_>::getLenght() const
{
    /**< 0 based coordinates, so N - N  is a legal fragment covering a single nucleotide at position N */
    return (m_endPos-m_startPos+1);
}

/** \brief Increase size of the element. Coordinates can go no lower then 0, but wiill not throw
 *
 * \param extendLeft int : Left shift size, must be +
 * \param extendRight int : Right shift size, must be +
 * \sa extendSite
 * \sa trimSite
 * \return void
 *
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::extendSite(long int extendLeft, long int extendRight)
{
    try
    {
        if((extendLeft<0)||(extendRight<0))
        {
            throw param_throw()<< string_error("INIT throwing from extendSite("
                                               +utility::to_string(extendLeft)+","+utility::to_string(extendRight)+"), param < 0 \n"  );;
        }
        int start=(getStart()-extendLeft);
        if (start < 0)
        {
            start=0;
        }
        setStart(start);
        setEnd(getEnd()+extendRight);
    }
    catch(param_throw & e )
    {
        elem_throw e;
        std::string * trace;


        if ( (trace=(boost::get_error_info<string_error>(e))) )
        {
            e << string_error(*trace+"Catching and re-throwing from extendSite("
                              +utility::to_string(extendLeft)+","+utility::to_string(extendRight)+")\n");
        }
        else
        {
            e << string_error("Catching and re-throwing from extendSite("
                              +utility::to_string(extendLeft)+","+utility::to_string(extendRight)+")\n");
            throw(e);
        }
    }
}

/** \brief Increase size of the element by the same for both begin and end. Coordinates can go no lower then 0,
 *
 * \param extend int : Size of shift on both sides
 * \return void
 *
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::extendSite(long int extend)
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
void uGenericNGS<_SELF_>::trimSite(long int trimLeft, long int trimRight)
{
    /**< Validate input */
    try
    {
        if ((trimLeft<0)||(trimRight<0)||(trimLeft+trimRight>this->getLenght()))
        {
            throw param_throw()<< string_error("PARAMERROR, throwing from trimSite("
                                               +utility::to_string(trimLeft)+","+utility::to_string(trimRight)+"), param < 0 \n"  );
        }
        this->m_startPos=(this->m_startPos+trimLeft);
        this->m_endPos=(this->m_endPos-trimRight);
    }
    catch (param_throw & err)
    {
#ifdef DEBUG
        std::cerr << "Throwing in trimSite(int, int)"<<std::endl;
#endif
        throw(err);
    }
}

/** \brief Diminish our element by an equal amount from both ends
 *
 * \param trim int : Size to trim from each side
 * \sa extendSite extendSite
 * \return void
 *
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::trimSite(long int trim)
{
    try
    {
        this->trimSite(trim,trim);
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
    {
        returnb=utility::isOverlap(this->m_startPos, this->m_endPos, other.m_startPos,other.m_endPos,type);
    }
    return returnb;
}

/** \brief Write contig to file using a writer.
 * \param uWriter& pWriter: A writer is used to specify the formating of the values and save them accordingly to disk.
 * \return void
 */
template <class _SELF_>
void uGenericNGS<_SELF_>::writeToOutput(uWriter& pWriter) const
{

    try
    {
        pWriter.writeToken(this->createToken());
    }
    catch(uWriter_exception_base &e)
    {
        addStringError(e,"Failed while calling writeToken in writeToOutput");
        throw e;
    }
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
        {
            throw param_throw() << string_error("Asking for more bins then lenght of Elem /n");
        }

        /**< If Strict and we cannot exactly fit our bins, fail */
        if ((ptype==SplitType::STRICT)&&(leftover!=0))
        {
            throw param_throw()<< string_error( "STRICT parameter not respected /n");
        }

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

            } // End of switch(ptype)
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
    catch(std::exception &e )
    {
#ifdef DEBUG
        std::cerr <<"Failed in divideIntoNBin()"<<std::endl;
#endif
        throw e ;
    }
    return returnVec;
}

/** \brief Divide the contig into contigs of size N
 * \param N int : Size of bins to create.
 * \param const SplitType: Determines how you managed leftover
 * \return vector<uGenericNGS> List of returned contigs
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
        {
            throw param_throw() << string_error("Asking for more bins then lenght of Elem /n");
        }

        /**< If Strict and we cannot exactly fit our bins, fail */
        if ((type==SplitType::STRICT)&&(leftover!=0))
        {
            throw param_throw()<< string_error( "STRICT parameter not respected /n");
        }

        /**< Once we pass the check, generate our basic bin */

        int curStart=getStart();
        for (int i=0; i<nbBin; i++)
        {
            returnVec.push_back(  _SELF_(getChr(),curStart, (curStart+(binSize-1))  ));
            curStart+=binSize;
        }

        /**< Then we deal with any leftovers according to policy */
        /**< If IGNORE, we do nothing */
        if(leftover!=0)
        {
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
            } // End of switch(type)
        }
    }
    catch(param_throw &e)
    {
        std::string ourtrace="Failed in divideIntoNBin() /n";
        if(  std::string const * trace=boost::get_error_info<string_error>(e) )
        {
            ourtrace+=*trace;
        }
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
    try
    {
        if (p_Pos>= ((int)m_score.size()))
        {
            m_score.resize(p_Pos+1);
        }
        m_score.at(p_Pos)=p_score;
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Fetch the score at a specific region,
 * \param int p_Pos: the position where the score should be fetched
 */
template <class _SELF_>
float uGenericNGS<_SELF_>::getScore(int p_Pos) const
{
    if (p_Pos>=((int)m_score.size()))
    {
        throw param_throw()<<string_error("Asked for Score position "+std::to_string(p_Pos)+"that is not set");
    }
    return m_score.at(p_Pos);
}

/** \brief Create the parser Token associated with Element
 *
 * \return uToken The token returned
 *
 */
template <class _SELF_>
uToken uGenericNGS<_SELF_>::createToken() const
{
    std::stringstream ss;
    ss << "CHR\t"<<this->getChr()<<"\nSTART_POS\t"<<this->getStart()<<"\n" << "END_POS\t"<<this->getEnd()<<"\n";

    if (getStrand()==StrandDir::FORWARD)
        ss << "STRAND\t"<<"+"<<"\n";
    else
        ss << "STRAND\t"<<"-"<<"\n";

    if (getScoreCount()>0)
    {
        ss << "SCORE\t"<<this->getScore()<<"\n";
    }
    try
    {
        return uToken(ss);
    }
    catch(uToken_exception_base &e)
    {
        addStringError(e, "Failed while creating token in uGenericNGS::getToken()");
        throw e;
    }
}

/**<  */

/** \brief Prints a human readable version of the element in no particular format.
 *
 * \param pOut std::ostream& Output to write to.
 * \return void
 *
 */
template <class _SELF_>
void  uGenericNGS<_SELF_>::print(std::ostream &pOut)const
{
    pOut<<"Chrom: "<<getChr()<<std::endl;
    pOut<<"Start: "<<utility::to_string(getStart())<<std::endl;
    pOut<<"End: " <<utility::to_string(getEnd())<<std::endl;

    if (m_score.size()>0)
    {
        pOut<<"Scores: ";
        for(auto value: m_score )
            pOut<<value<<" ";
        pOut<<std::endl;
    }

}

} // End of namespace NGS






/**< Legacy code, will return a uGenericNGS from  a tab delimited code. Deprecated, will be removed */
namespace factory
{
template<class _SELF_>
inline static _SELF_ makeNGSfromTabString(std::string tabString)
{
    utility::Tokenizer tabLine(tabString);

    std::string m_chrm;
    int start=0, end=0;
    try
    {
        tabLine.NextToken();
        m_chrm = tabLine.GetToken();
        tabLine.NextToken();
        start = utility::stoi(tabLine.GetToken());
        tabLine.NextToken();
        end = utility::stoi(tabLine.GetToken());

        return _SELF_(m_chrm,start,end);
    }
    catch(...)
    {
        throw;
    }
}
}




#endif // UFORMATBASE_H_INCLUDED
