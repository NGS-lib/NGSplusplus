#include <string>
#include <iostream>
#include <vector>
#include "uFormatBase.h"
#include <cassert>
#include <limits>
namespace NGS {
/** \brief Increase size of the element. Coordinates can go no lower then 0,
 *
 * \param extendLeft int : Left shift size, must be +
 * \param extendRight int : Right shift size, must be +
 * \return void
 *
 */
void uGenericNGS::extendSite(int extendLeft, int extendRight)
{
    try
    {
        if((extendLeft<0)||(extendRight<0))
            throw param_throw()<< string_error("INIT throwing from extendSite("+std::to_string(extendLeft)+","+std::to_string(extendRight)+"), param < 0 \n"  );;

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

    if(boost::get_error_info<generic_error>(e) ==NULL)
        e << generic_error(*this);
     if ( (trace=(boost::get_error_info<string_error>(e))) )
        e << string_error(*trace+"Catching and re-throwing from extendSite("+std::to_string(extendLeft)+","+std::to_string(extendRight)+")\n");
     else
         e << string_error("Catching and re-throwing from extendSite("+std::to_string(extendLeft)+","+std::to_string(extendRight)+")\n");
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
void uGenericNGS::extendSite(int extend)
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
void uGenericNGS::trimSites(int trimLeft, int trimRight)
{
    /**< Validate input */
    try
    {
        if ((trimLeft<0)||(trimRight<0)||(trimLeft+trimRight>this->getLenght()))
            throw param_throw()<< string_error("PARAMERROR, throwing from trimSites("+std::to_string(trimLeft)+","+std::to_string(trimRight)+"), param < 0 \n"  );

        this->startPos=(this->startPos+trimLeft);
        this->endPos=(this->endPos-trimRight);
    }
    catch (param_throw & err)
    {
        err << generic_error(*this);
        #ifdef DEBUG
           cerr << "Throwing in trimSites(int, int)"<<endl;
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
void uGenericNGS::trimSites(int trim)
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
bool uGenericNGS::doesOverlap(uGenericNGS other, OverlapType type) const
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
void uGenericNGS::writeBedToOuput(std::ostream &out, bool endLine) const
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
void uGenericNGS::writeToOutput(uWriter& pWriter) const
{
    pWriter.writeToken(this->getToken());
}



/** \brief Divide the contig into N contigs of equal size
 *
 * \param N int : Number of entries to create from the entry
 * \return vector<uGenericNGS> List of created contigs
 *
 */
std::vector<uGenericNGS> uGenericNGS::divideIntoNBin(int N,SplitType ptype)
{
    std::vector<uGenericNGS> returnVec;
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
            returnVec.push_back(  uGenericNGS(getChr(),curStart, (curStart+(binSize-1))  ));
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
        e << generic_error(*this);
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
std::vector<uGenericNGS> uGenericNGS::divideIntoBinofSize(const int N, const SplitType type)
{
    std::vector<uGenericNGS> returnVec;
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
        e << generic_error(*this);
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
};
/**< A score is an arbitray value set that can be used later */
/** \brief Set the score of a contig. Note that this involves resizing the vector, so settting arbitrarily large score counts can bust your memory
 * \param float p_score: the score to set
 * \param int p_Pos: the count of the score to set
 */
void uGenericNGS::setScore(float p_score, int p_Pos)
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
float uGenericNGS::getScore(int p_Pos) const
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
 uToken uGenericNGS::getToken() const
 {
    std::stringstream ss;
    ss << "CHR\t"<<this->getChr()<<"\nSTART_POS\t"<<this->getStart()<<"\n" << "END_POS\t"<<this->getEnd()<<"\n";
    if (getScore()!=std::numeric_limits<float>::infinity())
       ss << "SCORE\t"<<this->getScore()<<"\n";

    return uToken(ss);
 }


/**< Legacy code, will return a uGenericNGS from  a tab delimited code. Deprecated, will be removed */
namespace factory
{
uGenericNGS makeNGSfromTabString(std::string tabString)
{

    utility::Tokenizer tabLine(tabString);

    std::string chrm;
    int start=0, end=0;
    try
    {
        tabLine.NextToken();
        chrm = tabLine.GetToken();
        tabLine.NextToken();
        start = std::stoi(tabLine.GetToken());
        tabLine.NextToken();
        end = std::stoi(tabLine.GetToken());

        return uGenericNGS(chrm,start,end);
    }
    catch(...)
    {
        throw;
    }

}

}



} // End of namespace NGS
