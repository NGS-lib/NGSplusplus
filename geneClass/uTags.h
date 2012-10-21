#ifndef UTAGS_H_INCLUDED
#define UTAGS_H_INCLUDED

#include "uFormats.h"
#include <memory>
namespace NGS {
//Our Tag format
//We used this to store mapped tags from NGS experiments
//This is used for single End tags
class uTags: public uGenericNGS
{

#define FORWARCHARD '+'
#define REVERSECHAR '-'

private:
    //0 = FORWARD, 1 = REVERSE

    //Optional
    //Map score, quality of the alignement.
    short int mapScore=255;

    //TODO this is some heavy overhead considering we are very unlikely to user Phred, cigar and name.. store in a different configuration?
    //The actual sequence
    std::string sequence="";
    char* name=nullptr;
    char* phredScore=nullptr;
    char* cigar=nullptr;
    //Paired End
    // bool PE;
    //Unmapped
    bool Unmapped=true;
    //Using Samtools definition here. Replace above variables by this when possible.
    int flag=0;
    //Pointer to next element of the template, not sure we really want to use this?.
    uTags * pMate=nullptr;
    int PELenght=0;

public:

    uTags();
    uTags(uGenericNGS otherItem);
    uTags(uToken pToken);
    uTags(std::string pchr, int start, int end, StrandDir pstrand=StrandDir::FORWARD);
    uTags(const uTags& copy_from);
    uTags& operator=  (uTags const& assign_from);
    ~uTags();



    void writeBedToOuput(std::ostream &out) const;
    void writetoBedCompletePE( std::ostream &out);
    void writeCompletedPESamToOutput(std::ostream &out);
    bool writeTrimmedSamToOutput(std::ostream &out, int left, int right);
    void loadfromSamString(std::string peakInfo, bool minimal);
    void writeSamToOutput(std::ostream &out) const;

    void setMate(uTags* pTag);

    void debugElem() const;


    bool isMapped()
    {
        return Unmapped;
    }
    void setMapped(bool pmapped)
    {
        Unmapped=pmapped;
    }

    uTags* getMate()
    {
        return pMate;
    };

    void setCigar(std::string pcigar)
    {
        if (cigar!=nullptr)
        {
            delete []cigar;
            cigar = nullptr;
        }

        try
        {

            cigar = new char[pcigar.size()+1];

            int lenght=pcigar.copy(cigar,pcigar.size(),0 );
            cigar[lenght]='\0';

        }
        catch(std::exception alloc_except)
        {
            std::cerr <<"Failled allocation in setCigar()";
            throw;
        }
    };
    std::string getCigar() const
    {
        std::string returnStr;
        if (cigar==nullptr)
            returnStr="";
        else
            returnStr=cigar;
        return returnStr;
    };


    /** \brief Parse your flag and set the necessary values to stay coherent.
     *
     * \param pflag int
     * \return void
     *
     */
    void setFlag(int pflag)
    {
        flag=pflag;
    };
    int getFlag() const
    {
        return flag;
    };
    void setSequence(std::string pSeq)
    {
        sequence=pSeq;
    };

    std::string getSequence() const
    {
        return sequence;
    };

    void setPhred(std::string Phred)
    {
        if (phredScore!=nullptr)
        {
            delete []phredScore;
            phredScore = nullptr;
        }

        try
        {
            phredScore = new char[Phred.size()+1];

            int lenght=Phred.copy(phredScore,Phred.size(),0 );
            phredScore[lenght]='\0';
        }
        catch(std::exception alloc_except)
        {
            std::cerr <<"Failled allocation in setPhred()";
            throw;
        }


    };

    std::string getPhred() const
    {
        std::string returnStr;
        if (phredScore==nullptr)
            returnStr=="";
        else
            returnStr=phredScore;
        return returnStr;
    };


    void setName(std::string pName)
    {
        if (name!=nullptr)
        {
            delete []name;
            name = nullptr;
        }

        try
        {
            name= new char [pName.size()+1];
            int lenght = pName.copy(name, pName.size(),0);
            name[lenght]='\0';
        }
        catch(std::exception alloc_except)
        {
            std::cerr <<"Failled allocation in setName()";
            throw;
        }

    };

    std::string getName() const
    {
        std::string returnStr;
        if (name==nullptr)
            returnStr="";
        else
            returnStr=name;
        return returnStr;
    };

    bool isPE() const
    {
        //if lenght above 0, is paired end
        return PELenght;
    };

    void setPELenght(int lenght)
    {
        try
        {
            if  (lenght <0)
                throw 10;

            PELenght=lenght;
        }
        catch(...)
        {
            std::cerr <<"Negative Paired end Lenght"<<std::endl;
            throw;
        }
    }
    int getPeLenght() const
    {
        return PELenght;
    };

    void setMapQual(int score)
    {
        mapScore=score;
    }

    int getMapQual() const
    {
        return mapScore;
    };

};

//Generator specific to uTags (allow to set the name, temporary (hopefully) fix)
template <>
template <class T2>
uTags uGenericNGSChrom<uTags>::generateRandomSite
(const int size_,std::mt19937& engine,const uGenericNGSChrom<T2> &exclList, const int sigma, const std::string ID) const
{
    //TODO Sanity check here to make sure it is possible to generate the asked for tag.
    uTags returnTag;

    bool found=false;
    int size = size_;

    int max = this->getChromSize();

    while (!found)
    {
        {
            uTags temptag;
            if (sigma!=0)
            {
                std::normal_distribution<float> gaussian(size, sigma);
                size = (int)gaussian(engine);
            }

            if (size>=1)
            {
                int shift = size/2;
                //Generating our distribution at each call is probably needlesly heavy.. check to optimize this in time.
                std::uniform_int_distribution<int> unif((shift+1), (max-shift));
                int center = unif(engine);
                temptag.setEnd(center+shift);
                temptag.setStart(center-shift);
                if ((exclList.getSubset(temptag.getStart(),temptag.getEnd())).count()==0)
                {
                    found=true;
                    returnTag=temptag;
                    returnTag.setChr(this->getChr());
                    returnTag.setName(ID);

                }
            }
        }
    }

    return returnTag;
}

namespace factory
{
uTags makeTagfromSamString(std::string samString, bool minimal=false);
}


class uTagsChrom: public uGenericNGSChrom<uTags>
{

private:

public:

    uTagsChrom();
    uTagsChrom(std::string ourChrom);
    uTagsChrom(uGenericNGSChrom<uTags>);
    uTags getTag(int i)
    {
        return VecSites.at(i);
    };
    void matchPE();
    void writeTrimmedSamToOutput(std::ostream &out, int left, int right);
    void writetoBedCompletePE(std::ostream& out);
    void writeCompletedPESamToOutput(std::ostream &out);
    void writeSamToOutput(std::ostream &out) const;
    void writeSamHeaderLine(std::ostream &out) const;
    void outputBedFormat(std::ostream& out) const;

    std::vector<float> getRegionSignal(int start, int end, bool overlap);

    //TODO FIX THIS
 //   int findPrecedingSite(int position, int low, int high)
 //   {
      //  return uGenericNGSChrom<uTags>::findPrecedingSite(position, low, high);
 //   };

};

// TODO: Lot's of code that should move to parser?
/**< Our complete tag experiment */
/**< Taga can be very very big....how do we deal with this. */
/**< For now we don't, but at some point we will have to. */
class uTagsExperiment: public uGenericNGSExperiment<uTagsChrom, uTags>
{

private:
    //std::ifstream& samStream;

    // void loadSamStream(std::ifstream& ourStream){samStream =ourStream;};
    void parseSamHeader();
public:

    void loadFromSam(std::ifstream& ourStream, bool minimal= false);
    void loadFromSamWithParser(std::string);

    void loadSamHeader(std::ifstream& ourStream);
    void writeToBed(std::ostream& out) const;
    void setChromSize(std::string chrom, int size);
    void writetoBedCompletePE(std::ostream& out);
    void writeSamToOutput(std::ostream &out) const ;
    void writeCompletedPESamToOutput(std::ostream &out);
    void writeTrimmedSamToOutput(std::ostream &out, int left, int right);
    uTags nextSamLine(bool minimal=false);
    std::vector<float> getRegionSignal(std::string chrom, int start, int end, bool overlap);
};
} // End of namespace NGS
#endif // UTAGS_H_INCLUDED
