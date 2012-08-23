#ifndef UFORMATBASE_H_INCLUDED
#define UFORMATBASE_H_INCLUDED

#include <string>
#include <iostream>
#include <vector>
#include "utility/utility.h"
#include "uGeneException.h"
//Our basic Site for NGS format
//Very weak class as there are many differences in the functionality of derived classes.

enum class SplitType{STRICT, IGNORE, EXTEND, ADD};
class uGenericNGS
{

#define FORWARCHARD '+'
#define REVERSECHAR '-'

#define FORWARD 0
#define REVERSE 1


protected:
    std::string chr;
    long int startPos=0;
    long int endPos=0;
    bool strand=0;

public:

    /**< Constructor taking chromosome name, start and end */
    uGenericNGS(std::string pchr, int pstart, int pend):chr(pchr)
    {
        try
        {
            /**< Implicitly, our startPosition must always be before our end position */
            {
                setEnd(pend);
                setStart(pstart);
            }
        }
        catch( std::exception &e)
        {
            #ifdef DEBUG
            std::cerr << "Error in uGenericNGS(std::string pchr, int pstart, int pend). data is"<< pchr<< " "<< pend<<" "<< pstart << " "<<std::endl;
            #endif
            throw e;
        }
    };

    uGenericNGS()
    { };

    virtual ~uGenericNGS(){};
/**< End Constructor/Destructor */


    /**< Write to Bed functions */
    virtual void writeBedToOuput(std::ostream &out, bool endLine) const;
    virtual void writeBedToOuput(std::ostream &out) const {
        writeBedToOuput(out, true);
    };
    /**< Get /Set */
    std::string getChr() const
    {
        return chr;
    };
    void setChr(std::string ourchr)
    {
        chr=ourchr;
    };

    void setStart(int ourStart)
    {
        try
        {
            if (!((ourStart<=getEnd())&&(ourStart>=0)))
               throw 10;
            startPos=ourStart;

        }
        catch(int &x)
        {
            #ifdef DEBUG
                std::cerr << "throwing in setStart" <<std::endl;
            #endif
            elem_throw e;
            e << string_error("Failed in setStart, ourStart is smalled then end or under 0, start is "+utility::numberToString(ourStart)+ " end is "+ utility::numberToString((int)getEnd()) +"\n");
            throw e;
        }
    };
    void setEnd(int ourEnd)
    {
        try{
            if (!((ourEnd>=getStart())&&(ourEnd>=0)))
                throw 10;
            endPos=ourEnd;
        }
        catch(...)
        {
            #ifdef DEBUG
                std::cerr << "throwing in setEnd" <<std::endl;
            #endif
            elem_throw e;
            e << string_error("throwing in setEnd(), start at "+utility::numberToString((int)getStart())+ " end is "+ utility::numberToString(ourEnd) +"\n");
            throw e;
        }
    };
    //TODO make sure this offer appropriate guarantee
    void setStartEnd(long int ourStart, long int ourEnd){
        try{
            setEnd(ourEnd);
            setStart(ourStart);
        }
        catch(...)
        {
            std::cerr << "throwing in setStartEnd" <<std::endl;
            throw;
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
    #ifdef DEBUG
    virtual void debugElem() const
    {
        using namespace utility;
        stringTocerr("Outputting elemn data");
        stringTocerr("Chrom "+getChr());
        stringTocerr("Start "+numberToString((int)getStart()));
        stringTocerr("End " +numberToString((int)getEnd()));
    }
    #endif
    /**<  Divide our region into a certain number of subregions */

    std::vector<uGenericNGS> divideIntoBinofSize(const int N, const SplitType type=SplitType::STRICT);
    std::vector<uGenericNGS> divideIntoNBin(const int N, const SplitType ptype=SplitType::STRICT);

    //If we want to Change the dimensions of our site
    void extendSite(int extend);
    void extendSite(int extendLeft, int extendRight);

    void trimSites(int trim);
    void trimSites(int trimLeft,int trimRight);

    /**< Should this be there? */
    bool doesOverlap(uGenericNGS other,OverlapType type=OverlapType::OVERLAP_PARTIAL) const;


};

namespace factory{
  //  uTags makeTagfromSamString(std::string samString, bool minimal=false);
    uGenericNGS makeNGSfromTabString(std::string tabString);
}


#endif // UFORMATBASE_H_INCLUDED
