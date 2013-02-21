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
uGene::uFeature::uFeature(long long pStart, long long pEnd, featureType pType,std::string pID, std::string pClass, short int pOffset ): m_start(pStart), m_end(pEnd), m_type(pType)
    ,m_ID(pID),m_class(pClass),m_offset(pOffset)
{
    if ( pStart < 0 || pStart > pEnd )
        throw ugene_exception_base();
}

/**< uGene */

bool uGene::uFeature::operator==(const uFeature &other) const
{

    return ((this->m_start==other.m_start) &&
            (this->m_end == other.m_end) &&
            (this ->m_type == other.m_type) &&
            ( this->m_ID == other.m_ID) &&
            (this->m_class == other.m_class)&&
            (this->m_offset == other.m_offset));
}

bool uGene::uFeature::operator!=(const uFeature &other) const
{
    return !(*this == other);
}

/**< From uTokens */
uGene::uGene(uToken pToken)try:
    uGenericNGS(pToken)
{}
catch(ugene_exception_base &e)
{
#ifdef DEBUG
    std::cerr << "Error in uGene(uToken)." <<std::endl;
#endif
//TODO ugene error
    // e<<tag_error(*this);
    throw e;
}


uGene::uGene(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand):uGenericNGS(pChr,pStart,pEnd,pStrand)
{}

uGene::uGene(std::string pChr, long long int pStart, long long int pEnd, StrandDir pstrand, float pScore):uGenericNGS(pChr,pStart,pEnd,pstrand,pScore) {}
uGene::uGene(std::string pChr, long long int pStart, long long int pEnd, float pScore ):uGenericNGS(pChr,pStart,pEnd,pScore) {}



/** \brief Constructor taking a uTags as parameter.
 * \param uTags p_tags: the tag to add to create the uGene object from.
 */
uGene::uGene(uTags p_tags): uGenericNGS(p_tags.getChr(),p_tags.getStart(), p_tags.getEnd(), p_tags.getStrand())
{
    setScoreVector(p_tags.getScoreVector());
}

/** \brief Constructor taking a uRegion as parameter.
 * \param uRegion p_region: the region to add to create the uGene object from.
 */
uGene::uGene(uRegion p_region): uGenericNGS(p_region.getChr(),p_region.getStart(), p_region.getEnd(), p_region.getStrand())
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
//TODO
bool uGene::isEqual(const uGene & pCompared) const
{

    /**< Compare Features */
    if (this->m_featureVector.size()==pCompared.featureCount())
    {
        auto itrComp= pCompared.featureBegin();
        for(auto itr=this->m_featureVector.cbegin(); itr!=this->m_featureVector.cend(); itr++)
        {
            if (*itr != *itrComp)
                return false;
            itrComp++;
        }
    }
    else
        return false;

    return ((this->getChr()==pCompared.getChr())&&
            (this->getStrand()==pCompared.getStrand())&&
            (this->getStart()==pCompared.getStart())&&
            (this->getEnd()==pCompared.getEnd())&&
            (this->getScoreVector()==pCompared.getScoreVector())&&
            (this->getID()==pCompared.getID())&&
            (this->getClass()==pCompared.getClass())&&
            (this->getTranscript())==pCompared.getTranscript());
}

uGene uGene::getCopy()const
{
    uGene returnCopy = *this;
    return returnCopy;
}

void uGene::addFeature(long long pFeatureStart, long long pFeatureEnd, featureType pType,std::string pID, std::string pClass , short int pOffset)
{
    uFeature test(pFeatureStart,pFeatureEnd,pType, pID, pClass, pOffset);
    m_featureVector.push_back(test);
    stable_sort(m_featureVector.begin(),m_featureVector.end(),[](const uFeature & item1, const uFeature & item2)
    {
        return item1.getStart()<item2.getStart();
    });


    /**< Adjust boundary */
    if (pFeatureStart<m_BoundaryStart)
        m_BoundaryStart=pFeatureStart;
    if (pFeatureEnd > m_BoundaryEnd)
        m_BoundaryEnd = pFeatureEnd;
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

bool uGene::hasFeatureType(featureType pType)const
{
    return  ( std::find_if(m_featureVector.begin(), m_featureVector.end(), [&](const uFeature & item)
    {
        return item.getType()==pType;
    } )!=m_featureVector.end());
}

bool uGene::isOverlappingFeature(long long pStart, long long pEnd)
{
    for(auto & feature:m_featureVector)
    {
        if (utility::isOverlap(feature.getStart(),feature.getEnd(),pStart,pEnd))
            return true;
    }
    return false;
}
bool uGene::isOverlappingFeature(long long pStart, long long pEnd, featureType pType)
{
    for(auto & feature:m_featureVector)
    {
        if ( feature.getType()==pType && utility::isOverlap(feature.getStart(),feature.getEnd(),pStart,pEnd))
            return true;
    }
    return false;

}



featureType mapFeature(const std::string & pType )
{
    if (featureMap.count(pType))
        return featureMap.find(pType)->second;
    else
        return featureType::OTHER;

}

std::string featureStr(const featureType& pType)
{
    for (auto & curPar:featureMap)
    {
        if (curPar.second==pType)
            return curPar.first;

    }
    return "";
}


/**< uGeneChrom */


void uGeneChrom::addData(const uGene& newSite)
try
{
    if (newSite.getChr()!=this->chr)
        throw ugene_exception_base()<<string_error("adding base to Chrom with non-matching scaffold/chr name");
    /**< Check for ID/Transcript combination if an ID is set*/
    if (newSite.getID().size())
    {
        /**< If we are trying to add but the existing ID/Transcript combination is set, fail */
        if (this->getIDCount(newSite.getID(),newSite.getTranscript()));
        throw ugene_exception_base()<<string_error("adding base to uGeneChrom with already existing ID/Transcript pair");
    }
    m_isSorted=false;
    VecSites.push_back(newSite);
}
catch(ugene_exception_base & e)
{
    throw e;
}




/** \brief Parsing from a token is quite a bit more complex in this representation, as a given token may already be associated with a given ID or transcript.
 *
 * \param pToken const uToken
 * \return void
 *
 */
void uGeneChrom::addData(const uToken pToken)
{
    /**< We first validate if there is a ID and potentially a Transcript associated. If not, treat it as normal */
//std::cout <<pToken.getParam(token_param::CHR)<<"in "<<this->getChr()<<std::endl;
    if (pToken.getParam(token_param::CHR)!=this->getChr())
        throw ugene_exception_base()<<string_error("adding base to Chrom with non-matching scaffold/chr name");

    if (pToken.isParamSet(token_param::GROUP_ID))
    {
        std::string tokID="",tokTranscript="";
        /**<Validate if we have an existing region with the ID and Transcript.  */
        tokID=pToken.getParam(token_param::GROUP_ID);
        if (pToken.isParamSet(token_param::GROUP_TRANSCRIPT))
            tokTranscript= pToken.getParam(token_param::GROUP_TRANSCRIPT);
        std::vector<uGene>::iterator mainItr;
        int i=0;
        /**< Case where the ID is not  associated, we add first as main feature */
        if (getIDCount(tokID,tokTranscript)==0)
        {
            StrandDir dir=StrandDir::FORWARD;
            if ( pToken.isParamSet(token_param::STRAND) && pToken.getParam(token_param::STRAND)=="-")
                dir=StrandDir::REVERSE;

            uGene newItem(pToken.getParam(token_param::CHR),std::stoll(pToken.getParam(token_param::START_POS)),std::stoll(pToken.getParam(token_param::END_POS)),dir);
            if ( pToken.isParamSet(token_param::SCORE))
                newItem.setScore(std::stof(pToken.getParam(token_param::SCORE)));
            newItem.setID(tokID);
            newItem.setTranscript(tokTranscript);
            this->VecSites.push_back(newItem);
            this->m_isSorted=false;
            mainItr=VecSites.end();
            mainItr--;
            i++;
        }
        else
        {
            mainItr= this->findGene(tokID,tokTranscript);
        }
        for (; i<pToken.paramCount(token_param::START_POS); i++)
        {
            //TODO : Accessor
            StrandDir dir=StrandDir::FORWARD;
            // float score=0.0f;
            /**< Others as supplementary feature */

            if ( pToken.isParamSet(token_param::STRAND) && pToken.getParam(token_param::STRAND)=="-")
                dir=StrandDir::REVERSE;

            std::string featureID="", featureClass="", Myoffset;

            if (pToken.isParamSet(token_param::FEATURE_CLASS,i))
                featureClass=pToken.getParam(token_param::FEATURE_CLASS,i);
            if (pToken.isParamSet(token_param::FEATURE_ID))
                featureID=pToken.getParam(token_param::FEATURE_ID,i);
            if (pToken.isParamSet(token_param::OFFSET))
                featureID=pToken.getParam(token_param::FEATURE_ID,i);

            mainItr->addFeature(std::stoll(pToken.getParam(token_param::START_POS,i)),std::stoll(pToken.getParam(token_param::END_POS,i)),mapFeature(pToken.getParam(token_param::FEATURE_TYPE,i)),featureID,featureClass );
            if ( pToken.isParamSet(token_param::SCORE))
                mainItr->setScore(std::stof(pToken.getParam(token_param::SCORE)));

            mainItr->setID(tokID);
            mainItr->setTranscript(tokTranscript);

        }
    }
    else
    {
        /**< Otherwise normal group */
        uGenericNGSChrom::addData(pToken);
    }

}


long long uGeneChrom::getIDCount(const std::string & pId, const std::string & pTranscript)
{
    if (pTranscript.size()==0)
        return  ( std::count_if(std::begin(VecSites), std::end(VecSites), [&](const uGene & item)
    {
        return (item.getID()==pId) ;
    } ));
    else
        return  ( std::count_if(std::begin(VecSites), std::begin(VecSites), [&](const uGene & item)
    {
        return ( item.getID()==pId && item.getTranscript()==pTranscript) ;
    } ));
}

typename std::vector<uGene>::iterator uGeneChrom::findGene(const std::string & pId, const std::string & pTranscript)
{
    return std::find_if(std::begin(VecSites),std::end(VecSites),[&](uGene & item)
    {
        return(item.getID()==pId && item.getTranscript()== pTranscript ) ;
    });
}

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
    if (this == &copFrom)      // Same object?
        return *this;

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


/** \brief Similair to findNext, but will only find an item that contains a feature matching pType. Still works on current sort.
 *
 * \param pPosition long long  Position to from
 * \param pType featureType Feature type to filter on
 * \return typename std::vector<uGene>::const_iterator Iterator to the item, returns end() if not found.
 *
 */
typename std::vector<uGene>::const_iterator uGeneChrom::findNextWithFeature(long long pPosition, featureType pType)const
{
    try
    {
        /**< If unsorted, fail */

        if (VecSites.size()==0)
            return VecSites.end();

        if (m_isSorted==false)
            throw unsorted_throw() <<string_error("findNextGeneWithFeature called on unsorted vector \n") ;
        if ((sortGetStart==nullptr)||(sortGetEnd==nullptr))
            throw ugene_exception_base() <<string_error(" findNextGeneWithFeature called on chrom without appropriate start or end function\n") ;

        /**< Return true comparitor if value smaller then item 2 */
        auto comp = [&] (const float &posGiven, const uGene &item2)
        {

            return ( posGiven< sortGetStart(&item2) && item2.hasFeatureType(pType) );
        };
        /**< Compare, sort Value */
        auto upper = std::upper_bound(VecSites.begin(), VecSites.end(), pPosition, comp);

        /**< If no result*/
        if (upper==VecSites.end())
            return VecSites.end();

        /**<Found the upper bound, we now start iterating */
        return (upper);
    }
    catch (unsorted_throw & e)
    {
#ifdef DEBUG
        std::cerr << "findNextGeneWithFeature called on unsorted vector" <<std::endl;
#endif
        throw e;
    }
    catch (ugene_exception_base & e)
    {
#ifdef DEBUG
        std::cerr << "Calling findNextGeneWithFeature and you did not provide an aproriate get function" <<std::endl;
#endif
        throw e;
    }

}


/** \brief Similair to findPreceding, but will only find an item that contains a feature matching pType. Still works on current sort.
 *
 * \param pPosition long long  Position to from
 * \param pType featureType Feature type to filter on
 * \return typename std::vector<uGene>::const_iterator Iterator to the item, returns end() if not found.
 *
 */
typename std::vector<uGene>::const_iterator uGeneChrom::findPrecedingWithFeature(long long pPosition, featureType pType)const
{
    try
    {
        /**< If unsorted, fail */
        if (VecSites.size()==0)
            return VecSites.end();
        if (m_isSorted==false)
            throw unsorted_throw() <<string_error("findPrecedingSite called on unsorted vector \n") ;
        if ((sortGetStart==nullptr)||(sortGetEnd==nullptr))
            throw ugene_exception_base() <<string_error(" findPrecedingSite called on chrom without appropriate start or end function\n") ;
        auto comp = [&] (const uGene &item1, const float &pPos)
        {
            return ( sortGetStart(&item1)< pPos && item1.hasFeatureType(pType) );
        };

        /**< Compare, sort Value */
        auto lower = std::lower_bound(VecSites.begin(), VecSites.end(), pPosition, comp);
        /**< If result is our first item, then no item precedes it */
        if ((lower==VecSites.begin()))
            return VecSites.end();
        /**< If result is end, every idem precedes the value */
        /**<Return item precedes and as such is LESS then position. if no item was found, last item is closest to value  */
        lower--;
        return (lower);
    }
    catch (unsorted_throw & e)
    {
#ifdef DEBUG
        std::cerr << "FindPrecedingSite called on unsorted vector" <<std::endl;
#endif
        throw e;
    }
    catch (ugene_exception_base & e)
    {
#ifdef DEBUG
        std::cerr << "Calling findPrecedingSite and you did not provide an aproriate get function" <<std::endl;
#endif
        throw e;
    }
}

/** \brief Greedy search as our data is not sorted. In fact, we need to pass every element to check if it's boundary overlaps. Will search passed on current sort
 *
 * \param pPosition long long  Position to from
 * \param pType featureType Feature type to filter on
 * \return typename std::vector<uGene>::const_iterator Iterator to the item, returns end() if not found.
 *
 */
typename std::vector<uGene>::const_iterator uGeneChrom::findNextFeature(long long pPosition, featureType pType)const
{
    typename std::vector<uGene>::const_iterator curGene=VecSites.end();
    long long curPos=0;
    for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
    {
        if (it->hasFeatureType(pType) && it->getBoundaryEnd() > pPosition)
        {
            for (auto featureItr= it->featureBegin(); featureItr!=it->featureEnd(); featureItr++)
            {
                if (featureItr->getType()==pType && featureItr->getStart()>pPosition && featureItr->getStart()>curPos)
                {
                    curGene=it;
                    curPos= featureItr->getStart();
                }

            }
        }
    }
    return curGene;
}



/**< uGeneExperiment */
typename std::vector<uGene>::const_iterator uGeneExperiment::findNextGeneWithFeature(std::string pChr, long long pPosition, featureType pType)const
{
    try
    {
        if (!ExpMap.count(pChr))
        {
            throw param_throw()<<string_error("Failling in uGenericNGSExperiment::findPrecedingSite, value "+pChr+" does not exist.\n");
        }
        auto tempChrom = getpChrom(pChr);
        return tempChrom->findNextWithFeature(pPosition,pType);
    }
    catch(...)
    {
        throw;
    }

}

typename std::vector<uGene>::const_iterator uGeneExperiment::findPrecedingGeneWithFeature(std::string pChr,long long pPosition, featureType pType)const
{
    try
    {
        if (!ExpMap.count(pChr))
        {
            throw param_throw()<<string_error("Failling in uGenericNGSExperiment::findNextSite, value "+pChr+" does not exist.\n");
        }
        auto tempChrom = getpChrom(pChr);
        return tempChrom->findPrecedingWithFeature(pPosition,pType);
    }
    catch(...)
    {
        throw;
    }
}

//TODO validate this
void uGeneExperiment::addData(const uGeneChrom& pChrom)
{
    /**< If chrom Already exist,  */
    if (ExpMap.count(pChrom.getChr()) != 0)
    {
        uGeneChrom* currentChrom;
        currentChrom=&(ExpMap[pChrom.getChr()]);
        for (auto itChrom =pChrom.begin(); itChrom!= pChrom.end(); itChrom++)
        {
            currentChrom->addData(*itChrom);
        }
    }
    else /**< Make deep copy */
    {
        ExpMap.insert(std::pair<std::string,uGeneChrom>(pChrom.getChr(),pChrom));
    }
}

void uGeneExperiment::addData(uToken & pToken)
{

    try
    {
        std::string chr = pToken.getParam(token_param::CHR);

        uGeneChrom* ptempChrom;
        ptempChrom=&(ExpMap[chr]);
        ptempChrom->setChr(chr);
        ptempChrom->addData(pToken);
    }
    catch(std::exception & e)
    {
#ifdef DEBUG
        std::cerr << "Catching and re-throwing in uGeneExperiment::addData()" <<std::endl;
#endif
        throw ugene_exception_base();
    }

}


}
