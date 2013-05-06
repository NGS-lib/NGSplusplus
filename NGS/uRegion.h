// ***************************************************************************
// uRegion.h (c) 2013
// Alexei Nordell-Markovits : Sherbrooke University
// Charles Joly Beauparlant : Laval University
//
//       This file is part of the NGS++ library.
//
//    The NGS++ library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program (lgpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// ***************************************************************************


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
class uGene;
class uGeneChrom;
class uGeneExperiment;

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
    uRegion(uGene);

    ~uRegion();

    std::string getIdent() const
    {
        return ident;
    };

    void setIdent(std::string pIdent)
    {
        ident=pIdent;
    }

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

    void writeSignal(std::ostream& out, char pSep='\t')  const;
    bool isEqual(const uRegion & pCompared)const;

    uRegion getCopy()const;


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

class uRegionExperiment;
class uRegionChrom :  public uGenericNGSChrom<uRegionChrom,uRegion>
{
public:

    uRegionChrom():uGenericNGSChrom() {};
    uRegionChrom(std::string ourChr):uGenericNGSChrom(ourChr)
    { }
    uRegionChrom(std::string ourChr, long long int length):uGenericNGSChrom(ourChr,length)
    { }
    uRegionChrom(const uGenericNGSChrom<uRegionChrom,uRegion>&);
    uRegionChrom& operator=(const uRegionChrom& copFrom);
    uRegionChrom(const uRegionChrom&);
    uRegionChrom(const std::vector<uRegion> & copyVec):uGenericNGSChrom(copyVec){};

    uRegionChrom(uBasicNGSChrom);
    uRegionChrom(uTagsChrom);
    uRegionChrom(uGeneChrom);

    ~uRegionChrom(){};
    uRegionChrom getCopy()const;

    void measureDensityOverlap(const  uTagsChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const  uRegionChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const  uBasicNGSChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);
    void measureDensityOverlap(const  uGeneChrom& chromtoComp, const OverlapType=OverlapType::OVERLAP_PARTIAL);

    void generateSignal(const uGeneChrom & chromToComp);
    void generateSignal(const uRegionChrom & chromToComp);
    void generateSignal(const uTagsChrom & chromToComp);
    void generateSignal(const uBasicNGSChrom & chromToComp);

    void writeSignal(std::ostream& out, char pSep='\t');

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
    void measureDensityOverlap(const uGeneExperiment& expToComp, const OverlapType poverlap=OverlapType::OVERLAP_PARTIAL);

    void generateSignal(const uTagsExperiment& expToComp);
    void generateSignal(const uRegionExperiment & expToComp);
    void generateSignal(const uBasicNGSExperiment & expToComp);
    void generateSignal(const uGeneExperiment & expToComp);


    void writeSignal(std::ostream& out, char pSep='\t');

};
} // End of namespace NGS
#endif
