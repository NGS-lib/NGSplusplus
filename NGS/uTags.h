#ifndef UTAGS_H_INCLUDED
#define UTAGS_H_INCLUDED

#include "uFormats.h"
#include <memory>
#include "utility/utility.h"
namespace NGS {
//Our Tag format
//We used this to store mapped tags from NGS experiments
//This is used for single End tags
class uTags: public uGenericNGS<uTags>
{

#define FORWARCHARD '+'
#define REVERSECHAR '-'

private:
    //0 = FORWARD, 1 = REVERSE

    //Optional
    //Map score, quality of the alignement.
    short int mapScore=255; /*!<Total mapping quality of the read */

    std::string sequence=""; /*!<Sequence associated with the read. Optional */
    char* name=nullptr;  /*!<Name or ID associated with read. not guaranteed to be unique*/
    char* phredScore=nullptr; /*!<PhredScore associated with each position of the read*/
    char* cigar=nullptr; /*!<Cigar flag as defined by the SAM format*/
    bool Unmapped=true;
    //Using Samtools definition here. Replace above variables by this when possible.
    int flag=0;  /*!<Sam flag as defined by the SAM format*/
    //Pointer to next element of the template, not sure we really want to use this?.
    int PELenght=0; /**< Lenght of the total segment if paired. Can also be used for estimated lenght */

public:

    uTags();
    uTags(uGenericNGS otherItem);
    uTags(const uToken & pToken);
    uTags(std::string pChr, long long int pStart, long long int pEnd, float pScore);
    uTags(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD);
    uTags(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore);


    uTags(const uTags& copy_from);
    uTags& operator=  (uTags const& assign_from);
    ~uTags();


    void writetoBedCompletePE( std::ostream &out);
    void writeCompletedPESamToOutput(std::ostream &out);
    bool writeTrimmedSamToOutput(std::ostream &out, int left, int right);
    void loadfromSamString(std::string peakInfo, bool minimal);
    void writeSamToOutput(std::ostream &out) const;

    void setMate(uTags* pTag);

    void debugElem() const;


    /** \brief Returns true if the tag is unmapped or incorrectly mapped
     *
     * \return bool
     *
     */
    bool isMapped()
    {
        return Unmapped;
    }
    /** \brief Set the mapped value of the tag.
     *
     * \param pmapped bool : Value to set the status to
     * \return void
     *
     */
    void setMapped(bool pmapped)
    {
        Unmapped=pmapped;
    }


    /** \brief Set the cigar value
     *
     * \param pcigar std::string : String to set cigar to
     * \return void
     *
     */
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
        catch(std::exception & e)
        {
            #ifdef DEBUG
             std::cerr <<"Failled allocation in setCigar()";
            #endif
            throw  e;
        }
    };

    /** \brief get the Cigar string associated with the sequence
     *
     * \return std::string
     *
     */
    std::string getCigar() const
    {
        std::string returnStr;
        if (cigar==nullptr)
            returnStr="";
        else
            returnStr=cigar;
        return returnStr;
    };

    /** \brief Set the Sam Flag associated with the sequence
     *
     * \param pflag int : The sam flag to set
     * \return void
     *
     */
    void setFlag(int pflag)
    {
        flag=pflag;
    };
    /** \brief Get the Sam flag associated with the sequence.
     *
     * \return int : The Sam flag associated.
     *
     */
    int getFlag() const
    {
        return flag;
    };

    /** \brief Set the sequence associated with the element
     *
     *  Sets the sequence associated with the element. The sequence needs to be either
     *  null ("") or of size equal to the element.
     * \exception param_throw : Thrown when parameter size neighter null or equal to element getLenght()
     * \param pSeq std::string : The sequence to set
     * \return void
     *
     */
    void setSequence(std::string pSeq)
    {
        if(((int)pSeq.size()!=0)&&((int)pSeq.size()!=getLenght()))
            throw param_throw()<< string_error("Failling in seqSequence. Sequence size neither null or equal to element size.");
        sequence=pSeq;
    };
    /** \brief Get the sequence associated with the element.
     *
     * \return std::string : The sequence associated with the element.
     *
     */
    std::string getSequence() const
    {
        return sequence;
    };

    /** \brief Set the PhredScore assocaited with the sequence
     *
     * \param Phred std::string : The sequence to set
     * \exception
     * \return void
     *
     */
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
        catch(std::exception & e)
        {
            #ifdef DEBUG
            std::cerr <<"Failled allocation in setPhred()";
            #endif
            throw e;
        }


    };

    /** \brief get the PhredScore associated with the sequence
     *
     * \return std::string : The PhredScore vector
     *
     */
    std::string getPhred() const
    {
        std::string returnStr;
        if (phredScore==nullptr)
            returnStr=="";
        else
            returnStr=phredScore;
        return returnStr;
    };

    /** \brief Set the ID associated with the element.
     *
     *  This sets the ID associated with the element. There are no conditions attached to the ID
     *
     * \param pName std::string : ID to set the name to
     * \exception Thrown when new allocations fails.
     * \return void
     *
     */
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
        catch(std::exception & e)
        {
            #ifdef DEBUG
            std::cerr <<"Failled allocation in setName()";
            #endif

            throw e ;
        }

    };

    /** \brief Get the name/ID associated with the element
     *
     * \return std::string : The returned ID
     *
     */
    std::string getName() const
    {
        std::string returnStr;
        if (name==nullptr)
            returnStr="";
        else
            returnStr=name;
        return returnStr;
    };

    /** \brief True if the paired end lenght is above 0
     *
     * \return bool True if PElenght is set
     *
     */
    bool isPE() const
    {
        return PELenght;
    };

    /** \brief
     *
     * \param lenght int : Value to set PELenght to.
     * \exception : param_throw(): Throw if parameter is < 0.
     * \return void
     *
     */
    void setPELenght(int lenght)
    {

            if  (lenght <0)
                throw param_throw()<<string_error("Throwing in setPELenght. Set an invalid PE lenght<0");
            PELenght=lenght;

    }
    /** \brief Return the PELenght of the element
     *
     * \return int
     *
     */
    int getPeLenght() const
    {
        return PELenght;
    };

    /** \brief Set the mapping quality of the element, max value is 255
     *
     * \param score short int : Value to set MapQuality to
     * \return void
     *
     */
    void setMapQual(short int score)
    {
        mapScore=score;
    }

    /** \brief Return the mapping quality of the element.
     *
     * \return short int: Mapping value
     *
     */
    short int getMapQual() const
    {
        return mapScore;
    };

};

class uTagsChrom: public uGenericNGSChrom<uTagsChrom,uTags>
{

private:

public:

    uTagsChrom():uGenericNGSChrom(){};
    uTagsChrom(const std::string & ourChr):uGenericNGSChrom(ourChr)
    {
    }


    uTagsChrom(const uGenericNGSChrom<uTagsChrom,uTags>&);
    uTagsChrom& operator=(const uTagsChrom& copFrom);
    uTagsChrom(const uTagsChrom&);
    uTagsChrom(std::vector<uTags>);
    uTagsChrom(const std::vector<uTags> & copyVec):uGenericNGSChrom(copyVec){};
    uTags getTag(int i)
    {
        return VecSites.at(i);
    };
    template<class _OTHER_>
    uTags generateRandomSite(const int size_,std::mt19937& engine,const _OTHER_ &exclList, const int sigma, const std::string ID) const;
    void writeTrimmedSamToOutput(std::ostream &out, int left, int right);
    void writetoBedCompletePE(std::ostream& out);
    void writeCompletedPESamToOutput(std::ostream &out);
    void writeSamToOutput(std::ostream &out) const;
    void writeSamHeaderLine(std::ostream &out) const;
    void outputBedFormat(std::ostream& out) const;

    std::vector<float> getRegionSignal(int start, int end, bool overlap);

};

// TODO: Lot's of code that should move to parser?
/**< Our complete tag experiment */
class uTagsExperiment: public uGenericNGSExperiment<uTagsExperiment,uTagsChrom, uTags>
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

namespace factory
{
    NGS::uTags makeTagfromSamString(std::string samString, bool minimal=false);
}

#endif // UTAGS_H_INCLUDED
