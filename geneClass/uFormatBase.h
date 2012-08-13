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


//TODO : Replace member init with initialiser lists. More optimal
class uGenericNGS
{
protected:
    std::string chr;
    long int startPos=0;
    long int endPos=0;

public:
    uGenericNGS(std::string pchr, int pstart, int pend)
    {
        try
        {
            /**< Implicitly, our startPosition must always be before our end position */
            {
                setChr(pchr);
                setEnd(pend);
                setStart(pstart);
            }
        }
        catch(...)
        {
            std::cerr << "Error in uGenericNGS(std::string pchr, int pstart, int pend). data is"<< pchr<< " "<< pend<<" "<< pstart << " "<<std::endl;
            throw;
        }
    };

    uGenericNGS()
    {
        startPos=0;
        endPos=0;
    };

    virtual ~uGenericNGS(){};
/**< End Constructor/Destructor */

    virtual void writeBedToOuput(std::ostream &out) const {
        writeBedToOuput(out, true);
    };
    virtual void writeBedToOuput(std::ostream &out, bool endLine) const;
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
        catch(...)
        {
            std::cerr << "throwing in setStart" <<std::endl;
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
            std::cerr << "throwing in setEnd" <<std::endl;
            elem_throw e;
            e << string_error("throwing in setEnd(), start at "+utility::numberToString((int)getStart())+ " end is "+ utility::numberToString(ourEnd) +"\n");
            throw e;
        }
    };

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
    virtual void debugElem() const
    {
        using namespace utility;
        stringTocerr("Outputting elemn data");
        stringTocerr("Chrom "+getChr());
        stringTocerr("Start "+numberToString((int)getStart()));
        stringTocerr("End " +numberToString((int)getEnd()));
    }
    /**<  Divide our region into a certain number of subregions */
    std::vector<uGenericNGS> divideIntoBinofSize(int N, SplitType type=SplitType::STRICT);
    std::vector<uGenericNGS> divideIntoNBin(int N,SplitType ptype=SplitType::STRICT);

    //If we want to Change the dimensions of our site
    void extendSite(int extend);
    void extendSite(int extendLeft, int extendRight);

    void trimSites(int trim);
    void trimSites(int trimLeft,int trimRight);

    /**< Should this be there? */
    bool doesOverlap(uGenericNGS other,OverlapType type) const;
};


namespace factory{
  //  uTags makeTagfromSamString(std::string samString, bool minimal=false);
    uGenericNGS makeNGSfromTabString(std::string tabString);
}


#endif // UFORMATBASE_H_INCLUDED
