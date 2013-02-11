#include "uGene.h"
#include "uBasicNGS.h"
#include "uTags.h"
#include "uRegion.h"
namespace NGS
{



/** \brief Constructor for our internal feature representation
 *
 * \param pStart long long
 * \param pEnd long long
 * \param pType featureType
 * \param pOffset short int
 * \param pID std::string
 *
 */
uGene::uFeature::uFeature(long long pStart, long long pEnd, featureType pType, short int pOffset, std::string pID ): m_start(pStart), m_end(pEnd), m_type(pType)
    ,m_ID(pID),m_offset(pOffset) {}

/** \brief Constructor taking a uTags as parameter.
 * \param uTags p_tags: the tag to add to create the uGene object from.
 */
uGene::uGene(uTags p_tags)
    : uGenericNGS(p_tags.getChr(),p_tags.getStart(), p_tags.getEnd(), p_tags.getStrand())
{
    setScoreVector(p_tags.getScoreVector());
}

/** \brief Constructor taking a uRegion as parameter.
 * \param uRegion p_region: the region to add to create the uGene object from.
 */
uGene::uGene(uRegion p_region)
    : uGenericNGS(p_region.getChr(),p_region.getStart(), p_region.getEnd(), p_region.getStrand())
{
    setScoreVector(p_region.getScoreVector());
}
//TODO test this better
/** \brief Checks if the content of elements are identical
 *
 * \param pCompared const uGene: Item to compare
 * \return bool True if identical
 *
 */
bool uGene::isEqual(const uGene & pCompared) const
{

    return ((this->getChr()==pCompared.getChr())&&
            (this->getStrand()==pCompared.getStrand())&&
            (this->getStart()==pCompared.getStart())&&
            (this->getEnd()==pCompared.getEnd())&&
            (this->getScoreVector()==pCompared.getScoreVector()));
}

uGene uGene::getCopy()const
{
    uGene returnCopy = *this;
    return returnCopy;
}

void  uGene::addFeature(long long pFeatureStart, long long pFeatureEnd, featureType pType, short int pOffset, std::string pID)
{
    m_featureVector.push_back(uFeature(pFeatureStart,pFeatureEnd,pType, pOffset, pID));
    stable_sort(m_featureVector.begin(),m_featureVector.end(),[](const uFeature & item1, const uFeature & item2)
    {
        return item1.getStart()<item2.getStart();
    });

}
void uGene::removeFeature(std::vector<uFeature>::const_iterator pItrPos)
{
    try
    {
        m_featureVector.erase(utility::to_mutable_iterator(m_featureVector,pItrPos));
    }
    catch(...)
    {
        throw;
    }

}
void uGene::removeFeature(std::vector<uFeature>::const_iterator pStartItr,std::vector<uFeature>::const_iterator pEndItr)
{
    try
    {
        m_featureVector.erase(utility::to_mutable_iterator(m_featureVector,pStartItr),utility::to_mutable_iterator(m_featureVector,pEndItr));
    }
    catch(...)
    {
        throw;
    }
}

    typename std::vector<uGene::uFeature>::const_iterator uGene::featureBegin()const
    {
        return m_featureVector.cbegin();
    }
    typename std::vector<uGene::uFeature>::const_iterator uGene::featureEnd()const
    {
        return m_featureVector.cend();
    }

    bool uGene::hasFeatureType(featureType pType){
        return  ( std::find_if(m_featureVector.begin(), m_featureVector.end(), [&](const uFeature & item){return item.getType()==pType;} )!=m_featureVector.end());
    }


/**< uGeneChrom */
/** \brief Copy constructor
 * \param const uGeneChrom& initFrom: the uGeneChrom to copy.
 */
uGeneChrom::uGeneChrom(const uGeneChrom& initFrom)
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
 * \param const uGeneChrom& copFrom: the uGeneChrom to copy.
 */
uGeneChrom& uGeneChrom::operator=(const uGeneChrom& copFrom)
{
    VecSites=copFrom.returnVecData();
    chr= copFrom.getChr();
    m_isSorted=copFrom.m_isSorted;
    sortGetStart=copFrom.sortGetStart;
    sortGetEnd=copFrom.sortGetEnd;
    m_comptFunc=copFrom.m_comptFunc;
    chromSize=copFrom.chromSize;

    return *this;
}


uGeneChrom uGeneChrom::getCopy()const
{
    uGeneChrom returnCopy = *this;
    return returnCopy;
}


}
