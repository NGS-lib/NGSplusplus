#ifndef UBASICNGS_H_INCLUDED
#define UBASICNGS_H_INCLUDED

#include "uFormats.h"
#include <memory>
namespace NGS
{

class uTags;
class uRegion;
class uGene;
class uBasicNGS: public uGenericNGS<uBasicNGS>
{
public:


    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD);
    uBasicNGS():uGenericNGS(){}
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore );
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, float pScore );

    uBasicNGS(uToken pToken);
    uBasicNGS(uTags pTags);
    uBasicNGS(uRegion pRegion);
    uBasicNGS(uGene pGene);

    ~uBasicNGS(){};

    bool isEqual(const  uBasicNGS & pCompared) const;

    uBasicNGS getCopy()const;

};

class uBasicNGSChrom: public uGenericNGSChrom<uBasicNGSChrom,uBasicNGS>
{
public:
    uBasicNGSChrom():uGenericNGSChrom()
    {
    }

    uBasicNGSChrom(std::string chr):uGenericNGSChrom(chr)
    {
    }

    uBasicNGSChrom(std::string chr, long long int lenght):uGenericNGSChrom(chr,lenght)
    {
    }

  //  uBasicNGSChrom(const uGenericNGSChrom<uBasicNGSChrom,uBasicNGS> & copyCop);
    uBasicNGSChrom(const uBasicNGSChrom& initFrom);
    uBasicNGSChrom& operator=(const uBasicNGSChrom& copFrom);
    uBasicNGSChrom(const std::vector<uBasicNGS> & copyVec):uGenericNGSChrom(copyVec){};

    uBasicNGSChrom getCopy()const;
    ~uBasicNGSChrom(){};

};

class uBasicNGSExperiment : public uGenericNGSExperiment<uBasicNGSExperiment,uBasicNGSChrom,uBasicNGS>
{
    public:
    ~uBasicNGSExperiment(){};
    uBasicNGSExperiment& getCopy();
};

} // End of namespace NGS
#endif // UBASICNGS_H_INCLUDED
