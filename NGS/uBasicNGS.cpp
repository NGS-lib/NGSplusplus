#include "uBasicNGS.h"
#include "uTags.h"
#include "uRegion.h"
#include "uGene.h"
namespace NGS
{

/** \brief Constructor taking a uTags as parameter.
 * \param uTags p_tags: the tag to add to create the uBasicNGS object from.
 */
uBasicNGS::uBasicNGS(uTags pTags)
    : uGenericNGS(pTags.getChr(),pTags.getStart(), pTags.getEnd(), pTags.getStrand())
{
    setScoreVector(pTags.getScoreVector());
}

/** \brief Constructor taking a uRegion as parameter.
 * \param uRegion p_region: the region to add to create the uBasicNGS object from.
 */
uBasicNGS::uBasicNGS(uRegion pRegion)
    : uGenericNGS(pRegion.getChr(),pRegion.getStart(), pRegion.getEnd(), pRegion.getStrand())
{
    setScoreVector(pRegion.getScoreVector());
}

/** \brief Constructor taking a uRegion as parameter.
 * \param uRegion p_region: the region to add to create the uBasicNGS object from.
 */
uBasicNGS::uBasicNGS(uGene pGene)
    : uGenericNGS(pGene.getChr(),pGene.getStart(), pGene.getEnd(), pGene.getStrand())
{
    setScoreVector(pGene.getScoreVector());
}


/** \brief Constructor from parser Token
* \param uToken Valid Token with data. We assumed the token is valid, so skip some checks here to avoid duplication
*/
uBasicNGS::uBasicNGS(uToken pToken) try :
    uGenericNGS(pToken)
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
    uBasicNGS::uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand) try : uGenericNGS(pChr,pStart,pEnd,pStrand)
    { }
    catch(construct_elem_throw &e)
    {
        addStringError(e,"Throwing in uBasicNGS(string, long long int, long long int, Strandir)");
        e << basic_error(*this);
        throw e;
    }

    /** \brief Constructor with basic region infos plus strand and score
     * \param std::string chr: the name of the chromosome
     * \param long long int start: the starting position of the element to add.
     * \param long long int end: the ending position of the element to add.
     * \param StrandDir dir: the strand of the element to add.
     * \param float pScore: the score of the element to add.
     * \exception construct_elem_throw: When the start and/or end position are invalid.
     */
    uBasicNGS::uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore )
        : uGenericNGS(pChr, pStart, pEnd,pStrand,pScore)
    { }


    /** \brief Constructor with basic region infos plus score
     * \param std::string chr: the name of the chromosome
     * \param long long int start: the starting position of the element to add.
     * \param long long int end: the ending position of the element to add.
     * \param float pScore: the score of the element to add.
     * \exception construct_elem_throw: When the start and/or end position are invalid.
     */
    uBasicNGS::uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, float pScore )
        : uGenericNGS(pChr, pStart, pEnd,pScore)
    { }

/** \brief Checks if the content of elements are identical
 *
 * \param pCompared const &_SELF_ Item to compare
 * \return bool True if identical
 *
 */

bool uBasicNGS::isEqual(const uBasicNGS  &  pCompared) const
{

    return ((this->getChr()==pCompared.getChr())&&
            (this->getStrand()==pCompared.getStrand())&&
            (this->getStart()==pCompared.getStart())&&
            (this->getEnd()==pCompared.getEnd())&&
            (this->getScoreVector()==pCompared.getScoreVector()));
}


uBasicNGS uBasicNGS::getCopy()const
{
    uBasicNGS returnCopy = *this;
    return returnCopy;
}




/**< uBasicNGSChrom */

/** \brief Copy constructor
 * \param const uBasicNGSChrom& initFrom: the uBasicNGSChrom to copy.
 */
uBasicNGSChrom::uBasicNGSChrom(const uBasicNGSChrom& initFrom)
{
    VecSites=initFrom.returnVecData();
    chr= initFrom.getChr();
    m_isSorted=initFrom.m_isSorted;
    sortGetStart=initFrom.sortGetStart;
    sortGetEnd=initFrom.sortGetEnd;
    m_comptFunc=initFrom.m_comptFunc;
    chromSize=initFrom.chromSize;
}

/** \brief Assignment operator overload to copy a uBagicNGSChrom.
 * \param const uBasicNGSChrom& copFrom: the uBasicNGSChrom to copy.
 */
uBasicNGSChrom& uBasicNGSChrom::operator=(const uBasicNGSChrom& copFrom)
{

    if (this == &copFrom) return *this;
    VecSites=copFrom.returnVecData();
    chr= copFrom.getChr();
    m_isSorted=copFrom.m_isSorted;
    sortGetStart=copFrom.sortGetStart;
    sortGetEnd=copFrom.sortGetEnd;
    m_comptFunc=copFrom.m_comptFunc;
    chromSize=copFrom.chromSize;

    return *this;
}


uBasicNGSChrom uBasicNGSChrom::getCopy()const
{
    uBasicNGSChrom returnCopy = *this;
    return returnCopy;
}


}
