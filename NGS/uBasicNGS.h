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
public:
    /** \brief Constructor from parser Token
    * \param uToken Valid Token with data. We assumed the token is valid, so skip some checks here to avoid duplication
    */
uBasicNGS(uToken pToken) try : uGenericNGS(pToken)
    {}
    catch(construct_elem_throw &e)
    {
        addStringError(e,"Throwing in uBasicNGS(uToken)");
        e << basic_error(*this);
        throw e;
    }

    /** \brief Constructor overload with basic region infos with optional strand direction
     * \param std::string chr: the name of the chromosome
     * \param long long int start: the starting position of the element to add.
     * \param long long int end: the ending position of the element to add.
     * \param StrandDir dir: the strand of the element to add.
     * \exception construct_elem_throw: When the start and/or end position are invalid.
     */
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD) try : uGenericNGS(pChr,pStart,pEnd,pStrand)
    { }
    catch(construct_elem_throw &e)
    {
        addStringError(e,"Throwing in uBasicNGS(string, long long int, long long int, Strandir)");
        e << basic_error(*this);
        throw e;
    }

    /** \brief Empty constructor.
      */
    uBasicNGS():uGenericNGS()
    { }

    /** \brief Constructor with basic region infos plus strand and score
     * \param std::string chr: the name of the chromosome
     * \param long long int start: the starting position of the element to add.
     * \param long long int end: the ending position of the element to add.
     * \param StrandDir dir: the strand of the element to add.
     * \param float pScore: the score of the element to add.
     * \exception construct_elem_throw: When the start and/or end position are invalid.
     */
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore )
        : uGenericNGS(pChr, pStart, pEnd,pStrand,pScore)
    { }

    /** \brief Constructor with basic region infos plus score
     * \param std::string chr: the name of the chromosome
     * \param long long int start: the starting position of the element to add.
     * \param long long int end: the ending position of the element to add.
     * \param float pScore: the score of the element to add.
     * \exception construct_elem_throw: When the start and/or end position are invalid.
     */
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, float pScore )
        : uGenericNGS(pChr, pStart, pEnd,pScore)
    { }

    /** \brief Constructor using the parent class uGenericNGS
     * \param uGenericNGS otherItem: an object from the parent class of uBasicNGS, uGenericNGS.
     */
    uBasicNGS(uTags p_tags);
    uBasicNGS(uRegion p_region);
    bool isEqual(const  uBasicNGS & pCompared) const;

    uBasicNGS getCopy();

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


};

class uBasicNGSExperiment : public uGenericNGSExperiment<uBasicNGSExperiment,uBasicNGSChrom,uBasicNGS>
{

    uBasicNGSExperiment& getCopy();

};

} // End of namespace NGS
#endif // UBASICNGS_H_INCLUDED
