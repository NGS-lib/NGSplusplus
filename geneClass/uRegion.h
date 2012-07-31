#ifndef UREGION_H
#define UREGION_H

//A region, is a fairly generic entity
//General, we use this to measure statistics and overlaps on a given region..

#include "uFormats.h"
#include "uTags.h"


class uRegion : public uGenericNGS
{
    public:
        uRegion();
        uRegion(std::string chr, int start, int end);
        uRegion(uGenericNGS);
        virtual ~uRegion();

        std::string getIdent() const {return ident;};

        float getDensity() const {return density;};
        void setCount(int ucount){count=ucount;};
        int getCount()const {return count;};
        void setSignal(int i, float value);
        void setSignal(std::vector<float>);
        std::vector<float> getSignal();

        float getScore(int p_Pos) const;
        float getScore() const {return getScore(0);};
        int getScoreCount() const { return score.size();};


        void setScore(float p_score, int p_Pos);
        void setScore(float ourscore) {setScore(ourscore,0);}

        void writeSignal(std::ostream& out) const;
        void writeRegion(std::ostream& out) const;
        void writeAll(std::ostream& out ) const;
    protected:
    private:

    std::vector<float> signal;
    // User friendly name of our region.
    std::string ident;
    //Varies according to our needs.
    float density;
    //Varies again
    int count;
    //And again
    std::vector<float> score;
};

class uRegionExperiment;
class uRegionChrom :  public uGenericNGSChrom<uRegion>
{
    public:
     void measureDensityOverlap(uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>,uGenericNGS>& expToComp , OverlapType=OverlapType::OVERLAP_PARTIAL);
     void measureDensityOverlap(uGenericNGSChrom<uGenericNGS>& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
     void measureDensityOverlap(uTagsChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
     void writeDensityAsTab(std::ostream& out);
     void generateSignal(uTagsExperiment& expToComp);
     void writeAll(std::ostream& out);
     void writeSignal(std::ostream& out);
     void generateSignal(const uRegionExperiment & expToComp);
     void generateSignal(const uRegionChrom & chromToComp);

};


class uRegionExperiment: public uGenericNGSExperiment<uRegionChrom, uRegion>{
     public:
    void measureDensityOverlap(uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>, uGenericNGS>& expToComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(uTagsExperiment& expToComp, const OverlapType poverlap=OverlapType::OVERLAP_PARTIAL);
    void generateSignal(uTagsExperiment& expToComp);
    void generateSignal(const uRegionExperiment & expToComp);
    void writeDensityAsTab(std::ostream& out);
    void writeAll(std::ostream& out);
    void writeSignal(std::ostream& out);


   void loadFromWig(std::ifstream & inputStream);

};


#endif
