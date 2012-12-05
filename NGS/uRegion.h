#ifndef UREGION_H
#define UREGION_H

/**< A region, is a fairly generic entity */
/**< General, we use this to measure statistics and overlaps on a given region.. */

#include "uFormats.h"
#include "uTags.h"
#include <limits>
namespace NGS {
class uToken;
class uParser;
class uBasicNGS;
class uBasicNGSChrom;
class uBasicNGSExperiment;
class uRegion : public uGenericNGS
{
    public:
        uRegion();
        uRegion(std::string chr, int start, int end);
        uRegion(uGenericNGS);

        uRegion(uToken);

        virtual ~uRegion();

        std::string getIdent() const {return ident;};
		void setDensity(float pDensity) {density=pDensity ;};
        float getDensity() const {return density;};
        void setCount(int ucount) {count=ucount;};
        int getCount()const {return count;};
        void setSignal(int i, float value);
        void setSignal(std::vector<float>);
        std::vector<float> getSignal();


        void writeSignal(std::ostream& out) const;
        void writeRegion(std::ostream& out) const;
        void writeAll(std::ostream& out ) const;

        //TODO substract regions (with options)

    protected:
    private:

    std::vector<float> signal={};
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

     void measureDensityOverlap(const  uGenericNGSChrom<uGenericNGSChrom,uGenericNGS>& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
     void measureDensityOverlap(const  uTagsChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
     void measureDensityOverlap(const  uRegionChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
     void measureDensityOverlap(const  uBasicNGSChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);


     void generateSignal(const uRegionChrom & chromToComp);
     void generateSignal(const uTagsChrom & chromToComp);
     void generateSignal(const uBasicNGSChrom & chromToComp);

     void writeDensityAsTab(std::ostream& out);
     void writeAll(std::ostream& out);
     void writeSignal(std::ostream& out);

};

class uRegionExperiment: public uGenericNGSExperiment<uRegionExperiment,uRegionChrom, uRegion>{
     public:


    uRegionExperiment& operator=(const uRegionExperiment& copFrom)=default;
    uRegionExperiment(const uRegionExperiment&) = default;
    uRegionExperiment()=default;


    void measureDensityOverlap(const uTagsExperiment& expToComp, const OverlapType poverlap=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const uRegionExperiment& expToComp, const OverlapType poverlap=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const uBasicNGSExperiment& expToComp, const OverlapType poverlap=OverlapType::OVERLAP_PARTIAL);



    void generateSignal(const uTagsExperiment& expToComp);
    void generateSignal(const uRegionExperiment & expToComp);
    void generateSignal(const uBasicNGSExperiment & expToComp);

    void writeDensityAsTab(std::ostream& out);
    void writeAll(std::ostream& out);
    void writeSignal(std::ostream& out);

   void loadFromWig(std::ifstream & inputStream);
};
} // End of namespace NGS
#endif
