#ifndef UBASICNGS_H_INCLUDED
#define UBASICNGS_H_INCLUDED

#include "uFormats.h"
#include <memory>
namespace NGS
{

class uTags;
class uRegion;

class uBasicNGS: public uGenericNGS<uBasicNGS>
{
    /** \brief Constructor from parser Token
    *
    * \param uToken Valid Token with data. We assumed the token is valid, so skip some checks here to avoid duplication
    *
    */
public:
uBasicNGS(uToken pToken)try :
    uGenericNGS(pToken)
    {
    }
    catch(construct_elem_throw &e)
    {
        addStringError(e,"Throwing in uBasicNGS(uToken)");
        e << basic_error(*this);
        throw e;
    }

    uBasicNGS(std::string chr, long long int start, long long int end, StrandDir dir)try:uGenericNGS(chr,start,end,dir)
    {   }
    catch(construct_elem_throw &e)
    {
        addStringError(e,"Throwing in uBasicNGS(string, long long int, long long int, Strandir)");
        e << basic_error(*this);
        throw e;
    }
    uBasicNGS():uGenericNGS()
    { }

    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore ):uGenericNGS(pChr, pStart, pEnd,pStrand,pScore) {}
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, float pScore ):uGenericNGS(pChr, pStart, pEnd,pScore) {}
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD ):uGenericNGS(pChr, pStart, pEnd,pStrand) {}

    uBasicNGS(uGenericNGS otherItem):uGenericNGS(otherItem)
    { }
    uBasicNGS(uTags   p_tags);
    uBasicNGS(uRegion p_region);
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
    uBasicNGSChrom(const uGenericNGSChrom<uBasicNGSChrom,uBasicNGS> & copyCop);
    uBasicNGSChrom(const uBasicNGSChrom& initFrom);
    uBasicNGSChrom& operator=(const uBasicNGSChrom& copFrom);
    uBasicNGSChrom(const std::vector<uBasicNGS> & copyVec):uGenericNGSChrom(copyVec){};
};

class uBasicNGSExperiment: public uGenericNGSExperiment<uBasicNGSExperiment,uBasicNGSChrom,uBasicNGS>
{

};

}
#endif // UBASICNGS_H_INCLUDED
