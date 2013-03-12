#ifndef UREGION_H
#define UREGION_H

/**< A region, is a fairly generic entity */
#include "uFormats.h"
//#include "uTags.h"
#include <limits>
#include "utility/utility.h"
namespace NGS
{
class uToken;
class uParser;
class uBasicNGS;
class uBasicNGSChrom;
class uBasicNGSExperiment;
class uTags;
class uTagsChrom;
class uTagsExperiment;

class uRegion : public uGenericNGS<uRegion>
{
public:
    uRegion();
    uRegion(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD);
    uRegion(std::string pChr, long long int pStart, long long int pEnd, StrandDir pstrand, float pScore);
    uRegion(std::string pChr, long long int pStart, long long int pEnd, float pScore );

    uRegion(uTags);
    uRegion(uBasicNGS);
    uRegion(uToken);


    ~uRegion();

    std::string getIdent() const
    {
        return ident;
    };
    void setDensity(float pDensity)
    {
        density=pDensity ;
    };
    float getDensity() const
    {
        return density;
    };
    void setCount(int pCount)
    {
        count=pCount;
    };
    int getCount()const
    {
        return count;
    };
    void setSignal(int i, float value);
    void setSignal(std::vector<float>);
    std::vector<float> getSignal()const;


    void writeSignal(std::ostream& out) const;
 //  void writeRegion(std::ostream& out) const;
 //   void writeAll(std::ostream& out ) const;

    bool isEqual(const uRegion & pCompared)const;

    uRegion getCopy()const;

    //TODO substract regions (with options)

protected:
private:

    std::vector<float> signal= {};
    /**<  User friendly name of our region. */
    std::string ident="";
    /**< Varies according to our needs. */
    float density=std::numeric_limits<float>::infinity();
    /**< Varies again */
    int count=0;
    /**< And again */

};

// TODO: Move read-write to parser and output class?
class uRegionExperiment;
class uRegionChrom :  public uGenericNGSChrom<uRegionChrom,uRegion>
{
public:

    uRegionChrom():uGenericNGSChrom() {};
    uRegionChrom(std::string ourChr):uGenericNGSChrom(ourChr)
    { }
    uRegionChrom(std::string ourChr, long long int lenght):uGenericNGSChrom(ourChr,lenght)
    { }
    uRegionChrom(const uGenericNGSChrom<uRegionChrom,uRegion>&);
    uRegionChrom& operator=(const uRegionChrom& copFrom);
    uRegionChrom(const uRegionChrom&);
    uRegionChrom(const std::vector<uRegion> & copyVec):uGenericNGSChrom(copyVec){};

    uRegionChrom(uBasicNGSChrom);
    uRegionChrom(uTagsChrom);
    ~uRegionChrom(){};
    uRegionChrom getCopy()const;

    void measureDensityOverlap(const  uTagsChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const  uRegionChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const  uBasicNGSChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);

    void generateSignal(const uRegionChrom & chromToComp);
    void generateSignal(const uTagsChrom & chromToComp);
    void generateSignal(const uBasicNGSChrom & chromToComp);

//    void writeDensityAsTab(std::ostream& out);
 //   void writeAll(std::ostream& out);
    void writeSignal(std::ostream& out);

};

class uRegionExperiment: public uGenericNGSExperiment<uRegionExperiment,uRegionChrom, uRegion>
{
public:

    uRegionExperiment& operator=(const uRegionExperiment& copFrom)=default;
    uRegionExperiment(const uRegionExperiment&) = default;
    uRegionExperiment()=default;

    uRegionExperiment getCopy() const;
    ~uRegionExperiment(){};

    void measureDensityOverlap(const uTagsExperiment& expToComp, const OverlapType poverlap=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const uRegionExperiment& expToComp, const OverlapType poverlap=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const uBasicNGSExperiment& expToComp, const OverlapType poverlap=OverlapType::OVERLAP_PARTIAL);

    void generateSignal(const uTagsExperiment& expToComp);
    void generateSignal(const uRegionExperiment & expToComp);
    void generateSignal(const uBasicNGSExperiment & expToComp);

 //   void writeDensityAsTab(std::ostream& out);
//    void writeAll(std::ostream& out);
    void writeSignal(std::ostream& out);

  //  void loadFromWig(std::ifstream & inputStream);
};
} // End of namespace NGS
#endif
