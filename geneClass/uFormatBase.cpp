#include <string>
#include <iostream>
#include <vector>
#include "uFormatBase.h"
#include <cassert>
#include "utility.h"

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
            throw 20;

        int start=(getStart()-extendLeft);
        if (start < 0)
            start=0;

        setStart(start);
        setEnd(getEnd()+extendRight);

    }
    catch(elem_throw & e )
    {
        elem_throw e;
        std::string * trace;

    if(uGenericNGS const * errorTagPoint=boost::get_error_info<generic_error>(e) )
        e << generic_error(*this);
     if ( trace=(boost::get_error_info<string_error>(e) ))
        e << string_error(*trace+"Catching and re-throwing from extendSite("+utility::numberToString(extendLeft)+","+utility::numberToString(extendRight)+")\n");
        throw(e);
        return;
    }
    catch(int & err)
    {
        elem_throw e;
        e << string_error("INIT throwing from extendSite("+utility::numberToString(extendLeft)+","+utility::numberToString(extendRight)+"), param < 0 \n"  );
        e << generic_error(*this);
        throw(e);
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
    catch(int err)
    {
        throw;
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
            throw 20;
    }
    catch (int err)
    {
        std::cerr <<"Invalid parameters in trimSites(), nothing done, please validate"<<std::endl;
        throw(err);
    }

    this->startPos=(this->startPos+trimLeft);
    this->endPos=(this->endPos-trimRight);
}

/** \brief Diminish our element by an equal amount from both ends
 *
 * \param trim int : Size to trim from each side
 * \return void
 *
 */
void uGenericNGS::trimSites(int trim)
{
    this->trimSites(trim,trim);
}

/** \brief Verify if two sites overlap by at least 1 bp
 *
 * \param other uGenericNGS : Other site to validate
 * \return bool : True if overlap.
 *
 */
bool uGenericNGS::doesOverlap(uGenericNGS other) const
{
    bool returnb=false;
    if (getChr()==other.getChr())
        returnb=utility::checkOverlap(this->startPos, this->endPos, other.startPos,other.endPos);

    return returnb;
}


/** \brief Output in Bed format, Endline can be optionally be delegated
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

/** \brief Output in Bed format, Endline can be optionally be delegated
 *
 * \param out std::ostream& : Our output steram
 * \param endLine bool : If true, write endline
 * \return void
 *
 */

//void uGenericNGS::writeBedToOuput(std::ostream &out) const
//{
//    out << chr << "\t" << startPos << "\t" << endPos << std::endl;
//}



/** \brief
 *
 * \param N int : Number of entries to create from the entry
 * \return vector<uGenericNGS>
 *
 */
std::vector<uGenericNGS> uGenericNGS::divideIntoNBin(int N,SplitType ptype)
{
    std::vector<uGenericNGS> returnVec;
    int leftover = getLenght()%N;
    int binSize= getLenght()/N;
    try
    {



        if (this->getLenght()<N)
            throw(5);

        /**< If Strict and we cannot exactly fit our bins, fail */
        if ((ptype==SplitType::STRICT)&&(leftover!=0))
            throw 10;

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

            }
        }
    }
    catch(int err)
    {
        switch (err){

        case 5:
             std::cerr <<"Failed on  more bins then possible"<<std::endl;
            std::cerr <<"Lenght is =  "<< this->getLenght() << "Bins = " << N  << std::endl;
                break;
        case 10:
             std::cerr <<"Failed on STRICT and leftover >0"<<std::endl;
             std::cerr <<"Leftover =  "<< leftover << "Type = " << (int)ptype  << std::endl;
                break;
        }

        std::cerr <<"Failed in divideIntoNBin()"<<std::endl;
        throw;
    }
    catch( ... ){
        std::cerr <<"Failed in divideIntoNBin()"<<std::endl;
        throw;
    }

    return returnVec;
}

/** \brief
 *
 * \param N int : Size of bins to create.
 * \return vector<uGenericNGS>
 *
 */
std::vector<uGenericNGS> uGenericNGS::divideIntoBinofSize(int N, SplitType type)
{
    std::vector<uGenericNGS> returnVec;
    try
    {
        int binSize= N;
        int nbBin= getLenght()/N;
        int leftover = getLenght()%binSize;

        if (getLenght()<N)
            throw(10);

        /**< If Strict and we cannot exactly fit our bins, fail */
        if ((type==SplitType::STRICT)&&(leftover!=0))
            throw 10;

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

            }
    }
    catch(...)
    {
        std::cerr <<"Failed in divideIntoNBin()"<<std::endl;
        throw;
    }

    return returnVec;
};


namespace factory
{
/** @brief Create a basic format from a tab delimited formatted string
  *
  * @todo: document this function
  */
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
        start = utility::stringToInt(tabLine.GetToken());
        tabLine.NextToken();
        end = utility::stringToInt(tabLine.GetToken());

        return uGenericNGS(chrm,start,end);
    }
    catch(...)
    {
        throw;
    }

}

}

