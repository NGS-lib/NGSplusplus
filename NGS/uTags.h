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
    char* m_cigar=nullptr; /*!<Cigar flag as defined by the SAM format*/
    bool Unmapped=true;
    //Using Samtools definition here. Replace above variables by this when possible.
    int flag=0;  /*!<Sam flag as defined by the SAM format*/
    //Pointer to next element of the template, not sure we really want to use this?.
    int PELenght=0; /**< Lenght of the total segment if paired. Can also be used for estimated lenght */

public:

    uTags();
    uTags(const uBasicNGS & otherItem);
    uTags(const uRegion & otherItem);
    uTags(uToken pToken);

    uTags(std::string pChr, long long int pStart, long long int pEnd, float pScore);
    uTags(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD);
    uTags(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore);

    uTags(const uTags& copy_from);
    uTags& operator=  (uTags const& assign_from);
    ~uTags();

    uTags getCopy()const;
    uTags getCompletedCopy()const;


    bool isEqual(const uTags & pCompared)const;

    //TODO

    uToken createToken()const;


 //   void writetoBedCompletePE( std::ostream &out);
 //   void writeCompletedPESamToOutput(std::ostream &out);
 //   bool writeTrimmedSamToOutput(std::ostream &out, int left, int right);
 //   void loadfromSamString(std::string peakInfo, bool minimal);
 //   void writeSamToOutput(std::ostream &out) const;

    void print(std::ostream &pOut) const; /**< Prints a human readable version of the data */

    bool isMapped() const;
    void setMapped(bool pmapped);

    void setCigar(std::string pCigar);
    std::string getCigar() const;

    void setFlag(int pflag);
    int getFlag() const;

    void setSequence(std::string pSeq);
    std::string getSequence() const;

    void setPhred(std::string Phred);
    std::string getPhred() const;

    void setName(std::string pName);

    std::string getName() const;
    bool isPE() const;

    void setPELenght(int lenght);
    int getPeLenght() const;

    void setMapQual(short int score);

    short int getMapQual() const;
};


class uRegionChrom;
class uBasicNGSChrom;
class uTagsChrom: public uGenericNGSChrom<uTagsChrom,uTags>
{

private:

public:

    uTagsChrom():uGenericNGSChrom(){};
    uTagsChrom(const std::string & ourChr):uGenericNGSChrom(ourChr)
    { }

    uTagsChrom getCopy() const;

    uTagsChrom(const uGenericNGSChrom<uTagsChrom,uTags>&);
    uTagsChrom& operator=(const uTagsChrom& copFrom);
    uTagsChrom(const uTagsChrom&);
    uTagsChrom(const uRegionChrom &);
    uTagsChrom(const uBasicNGSChrom &);


    uTagsChrom(const std::vector<uTags> & copyVec):uGenericNGSChrom(copyVec){};
   // uTags getTag(int i)
  //  {
  //      return VecSites.at(i);
  //  };
    template<class _OTHER_>
    uTags generateRandomSite(const int size_,std::mt19937& engine,const _OTHER_ &exclList, const int sigma, const std::string ID) const;
 //   void writeTrimmedSamToOutput(std::ostream &out, int left, int right);
 //   void writetoBedCompletePE(std::ostream& out);
 //   void writeCompletedPESamToOutput(std::ostream &out);
 //   void writeSamToOutput(std::ostream &out) const;
 //   void writeSamHeaderLine(std::ostream &out) const;
 //   void outputBedFormat(std::ostream& out) const;

    std::vector<float> getRegionSignal(int start, int end, bool overlap);

};

// TODO: Lot's of code that should move to parser?
/**< Our complete tag experiment */
class uTagsExperiment: public uGenericNGSExperiment<uTagsExperiment,uTagsChrom, uTags>
{

private:
    //std::ifstream& samStream;

    // void loadSamStream(std::ifstream& ourStream){samStream =ourStream;};
  //  void parseSamHeader();
public:
   // void loadFromSam(std::ifstream& ourStream, bool minimal= false);
 //   void loadFromSamWithParser(std::string);
 //   void loadSamHeader(std::ifstream& ourStream);
 //   void writeToBed(std::ostream& out) const;
 //   void setChromSize(std::string chrom, int size);
 //   void writetoBedCompletePE(std::ostream& out);
 //   void writeSamToOutput(std::ostream &out) const ;
 //   void writeCompletedPESamToOutput(std::ostream &out);
 //   void writeTrimmedSamToOutput(std::ostream &out, int left, int right);
 //   uTags nextSamLine(bool minimal=false);
    std::vector<float> getRegionSignal(std::string chrom, int start, int end, bool overlap);

    uTagsExperiment getCopy() const;

};
} // End of namespace NGS

namespace factory
{
    NGS::uTags makeTagfromSamString(std::string samString, bool minimal=false);
}

#endif // UTAGS_H_INCLUDED
