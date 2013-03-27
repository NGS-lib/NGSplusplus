#include "uGene.h"
#include "uBasicNGS.h"
#include "uTags.h"
#include "uRegion.h"
namespace NGS
{

/** \brief Constructor for our internal feature representation
 *
 * \param pStart long int
 * \param pEnd long int
 * \param pEnd long int
 * \param pType featureType
 * \param pOffset short int
 * \param pID std::stringF
 *
 */
uFeature::uFeature(long int pStart, long int pEnd,StrandDir pDir, featureType pType,std::string pID, short int pOffset ): m_start(pStart), m_end(pEnd), m_type(pType)
    ,m_ID(pID),m_offset(pOffset)
{
    if ( pStart < 0 || pStart > pEnd )
        throw ugene_exception_base();
    setStrand(pDir);
}

uToken uGene::createToken() const
{
    std::stringstream ss;
    ss << "CHR\t"<<this->getChr()<<"\nSTART_POS\t"<<this->getStart()<<"\n" << "END_POS\t"<<this->getEnd()<<"\n";
    ss<<"FEATURE_TYPE\t"<<"GENE\n";
    if (getStrand()==StrandDir::FORWARD)
        ss << "STRAND\t"<<"+"<<"\n";
    else
        ss << "STRAND\t"<<"-"<<"\n";

    if (getScoreCount()>0)
    {
        ss << "SCORE\t"<<this->getScore()<<"\n";
    }
    /**< Gene does not offset */
     ss << "PHASE\t"<<"0"<<"\n";

    if (featureCount()>0)
    {
        for(auto itr=featureBegin(); itr!=featureEnd(); itr++)
        {
            ss << "START_POS\t"<<itr->getStart()<<"\n";
            ss << "END_POS\t"<<itr->getEnd()<<"\n";
            ss << "FEATURE_TYPE\t"<<featureString(itr->getType())<<"\n";
            ss << "PHASE\t"<<itr->getOffSet()<<"\n";
            if (itr->getID()!="")
            {
                ss << "FEATURE_ID\t"<<(itr->getID())<<"\n";
            }
        }

    }


    try
    {
        return uToken(ss);
    }
    catch(uToken_exception_base &e)
    {
        addStringError(e, "Failed while creating token in uGene::getToken()");
        throw e;
    }


}

/**< uGene */

/** \brief Equality operator overloadeto compare features
 *
 * \param other const uFeature&:  Feature to compare with
 * \return bool Comparison result
 *
 */
bool uFeature::operator==(const uFeature &other) const
{

    return ((this->m_start==other.m_start) &&
            (this->m_end == other.m_end) &&
            (this ->m_type == other.m_type) &&
            ( this->m_ID == other.m_ID) &&
            //  (this->m_class == other.m_class)&&
            (this->m_offset == other.m_offset));
}

/** \brief Not equal operator for Features
 *
 * \param other const uFeature& Feature to compare with
 * \return bool Result of the comparison
 *
 */
bool uFeature::operator!=(const uFeature &other) const
{
    return !(*this == other);
}


/** \brief Constructor taking a uToken. If multiple elements/features are present, the constructor will try to parse and add the additional features
 *
 * \param pToken uToken
 *
 */
uGene::uGene(uToken pToken)
{
try {
        std::string tokID="",tokTranscript="",myOffset="0";

        /**<Validate if we have an existing region with the ID and Transcript.  */
        if (pToken.isParamSet(token_param::GROUP_ID))
            tokID=pToken.getParam(token_param::GROUP_ID);
        if (pToken.isParamSet(token_param::GROUP_TRANSCRIPT))
            tokTranscript= pToken.getParam(token_param::GROUP_TRANSCRIPT);

        /**< Case where the ID is not  associated, we add first as main feature */
        StrandDir dir=StrandDir::FORWARD;
        if ( pToken.isParamSet(token_param::STRAND) && pToken.getParam(token_param::STRAND)=="-")
            dir=StrandDir::REVERSE;
        this->setStrand(dir);
        this->setChr(pToken.getParam(token_param::CHR));
        this->setStartEnd(std::stoll(pToken.getParam(token_param::START_POS)),std::stoll(pToken.getParam(token_param::END_POS)));

        if ( pToken.isParamSet(token_param::SCORE))
            this->setScore(std::stof(pToken.getParam(token_param::SCORE)));
        this->setID(tokID);
        this->setTranscript(tokTranscript);

        for (int i=1; i<pToken.paramCount(token_param::START_POS); i++)
        {
            std::string featureID;
            myOffset="0";
            StrandDir dir=StrandDir::FORWARD;
            /**< Others as supplementary feature */
            if ( pToken.isParamSet(token_param::STRAND,i) && pToken.getParam(token_param::STRAND,i)=="-"){
                dir=StrandDir::REVERSE;}
            if (pToken.isParamSet(token_param::FEATURE_ID,i)){
                featureID=pToken.getParam(token_param::FEATURE_ID,i); }
            if (pToken.isParamSet(token_param::PHASE,i)){
                myOffset=pToken.getParam(token_param::PHASE,i);}

            featureType curType=featureType::OTHER;
            if (pToken.isParamSet(token_param::FEATURE_TYPE,i)){
                    curType=mapFeature(pToken.getParam(token_param::FEATURE_TYPE,i));
                    }

            this->addFeature(std::stoll(pToken.getParam(token_param::START_POS,i)),std::stoll(pToken.getParam(token_param::END_POS,i)),dir,curType,featureID,std::stoi(myOffset));

        }

}
catch(...)
{
#ifdef DEBUG
    std::cerr << "Error in uGene(uToken)." <<std::endl;
#endif
    throw;
}
}

uGene::uGene(std::string pChr, long int pStart, long int pEnd, StrandDir pStrand):uGenericNGS(pChr,pStart,pEnd,pStrand)
{}

uGene::uGene(std::string pChr, long int pStart, long int pEnd, StrandDir pstrand, float pScore):uGenericNGS(pChr,pStart,pEnd,pstrand,pScore) {}
uGene::uGene(std::string pChr, long int pStart, long int pEnd, float pScore ):uGenericNGS(pChr,pStart,pEnd,pScore) {}

uGene::uGene(std::string pChr, long int pStart, long int pEnd, StrandDir pStrand, std::string pID, std::string pTranscript):uGenericNGS(pChr,pStart,pEnd,pStrand),m_ID(pID),m_transcript(pTranscript)
{}


    uGene::uGene(std::string pChr, long int pStart, long int pEnd, StrandDir pStrand, float pScore, std::string pID, std::string pTranscript):uGenericNGS(pChr,pStart,pEnd,pStrand,pScore),m_ID(pID),m_transcript(pTranscript)
    {
    }

    uGene::uGene(std::string pChr, long int pStart, long int pEnd, std::string pID, std::string pTranscript):uGenericNGS(pChr,pStart,pEnd),m_ID(pID),m_transcript(pTranscript)
    {
    }




/** \brief Constructor taking a uBasicNGS as parameter.
 * \param uTags p_tags: the tag to add to create the uGene object from.
 */

uGene::uGene(uBasicNGS pBasic): uGenericNGS(pBasic.getChr(),pBasic.getStart(), pBasic.getEnd(), pBasic.getStrand())
{
    setScoreVector(pBasic.getScoreVector());
}

/** \brief Constructor taking a uTags as parameter.
 * \param uTags p_tags: the tag to add to create the uGene object from.
 */
uGene::uGene(uTags pTag): uGenericNGS(pTag.getChr(),pTag.getStart(), pTag.getEnd(), pTag.getStrand())
{
    setScoreVector(pTag.getScoreVector());
}

/** \brief Constructor taking a uRegion as parameter.
 * \param uRegion p_region: the region to add to create the uGene object from.
 */
uGene::uGene(uRegion pRegion): uGenericNGS(pRegion.getChr(),pRegion.getStart(), pRegion.getEnd(), pRegion.getStrand())
{
    setScoreVector(pRegion.getScoreVector());
}

/** \brief Checks if the content of elements are identical
 *
 * \param pCompared const uGene: Item to compare
 * \return bool True if identical
 *
 */
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
    /**< Compare every other element */
    return ((this->getChr()==pCompared.getChr())&&
            (this->getStrand()==pCompared.getStrand())&&
            (this->getStart()==pCompared.getStart())&&
            (this->getEnd()==pCompared.getEnd())&&
            (this->getScoreVector()==pCompared.getScoreVector())&&
            (this->getID()==pCompared.getID())&&
            //  (this->getClass()==pCompared.getClass())&&
            (this->getTranscript())==pCompared.getTranscript());
}

/** \brief Return a deep copy of this element.
 *
 * \return uGene Copy of current element
 *
 */
uGene uGene::getCopy()const
{
    uGene returnCopy = *this;
    return returnCopy;
}

/** \brief Add a specified feature and associated with the uGene element. Also sets boundary as a byproduct if needed
 *
 * \param pFeatureStart long int: Start position of the feature, must be > then End pos
 * \param pFeatureEnd long int: End position of the feature
 * \param pStrand StrandDir : Direction of the feature
 * \param pType featureType : Feature type from valid list
 * \param pID std::string : Feature ID
 * \param pClass std::string : Feature Class
 * \param pOffset short int : Offset, typically set to 0-1-2
 * \return void
 *
 */
void uGene::addFeature(long int pFeatureStart, long int pFeatureEnd,StrandDir pStrand, featureType pType,std::string pID , short int pOffset)
{
    m_featureVector.push_back(uFeature(pFeatureStart,pFeatureEnd,pStrand,pType, pID, pOffset));
    /**< Resort features */
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


/** \brief Private utility function, sets the boundary to the min/max of the features
 *
 * \return void
 *
 */
void uGene::inferBoundary()
{
    m_BoundaryStart=(std::min_element(m_featureVector.begin(),m_featureVector.end(), [](const uFeature & isSmaller, const uFeature & isBigger)
    {
        return (isSmaller.getStart()<isBigger.getStart()) ;
    }))->getStart();
    m_BoundaryEnd=(std::max_element(m_featureVector.begin(),m_featureVector.end(),[](const uFeature & isSmaller, const uFeature & isBigger)
    {
        return (isSmaller.getEnd()<isBigger.getEnd()) ;
    }))->getEnd();
}

/** \brief Erase the feature pointed to by the iterator and resize boundary as needed
 *
 * \param pItrPos std::vector<uFeature>::const_iterator Item pointed to
 * \return void
 *
 */
void uGene::removeFeature(std::vector<uFeature>::const_iterator pItrPos)
{
    try
    {
        m_featureVector.erase(utility::to_mutable_iterator(m_featureVector,pItrPos));
        this->inferBoundary();
    }
    catch(...)
    {
        throw;
    }

}
/** \brief Erase range [) pointed to by two iterators, resize boundary as needed
 *
 * \param pStartItr std::vector<uFeature>::const_iterator Begining of range
 * \param pEndItr std::vector<uFeature>::const_iterator End of range, will not be erased.
 * \return void
 *
 */
void uGene::removeFeature(std::vector<uFeature>::const_iterator pStartItr,std::vector<uFeature>::const_iterator pEndItr)
{
    try
    {
        m_featureVector.erase(utility::to_mutable_iterator(m_featureVector,pStartItr),utility::to_mutable_iterator(m_featureVector,pEndItr));
        this->inferBoundary();
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Return const iterator to the first feature
 *
 * \return typename std::vector<uGene::uFeature>::const_iterator
 *
 */
typename std::vector<uFeature>::const_iterator uGene::featureBegin()const
{
    return m_featureVector.cbegin();
}

/** \brief Return const iterator to end of feature vector
 *
 * \return typename std::vector<uGene::uFeature>::const_iterator
 *
 */
typename std::vector<uFeature>::const_iterator uGene::featureEnd()const
{
    return m_featureVector.cend();
}

/** \brief Verify if the uGene has a feature with the given featureType
 *
 * \param pType featureType: Type to verify
 * \return bool: True if featureType found
 *
 */
bool uGene::hasFeatureType(featureType pType)const
{

    auto itemEqual =[&](const uFeature & item)
    {
        return item.getType()==pType;
    };

    return  ( std::find_if(m_featureVector.begin(), m_featureVector.end(), itemEqual ) !=m_featureVector.end() );


}

/** \brief Verifies if a feature specifically overlaps the given range. Only valide for positions
 *
 * \param pStart long int: Begining of range
 * \param pEnd long int : End of range
 * \return bool: True if at least one overlap found
 *
 */
bool uGene::isOverlappingFeature(long int pStart, long int pEnd)
{
    for(auto & feature:m_featureVector)
    {
        if (utility::isOverlap(feature.getStart(),feature.getEnd(),pStart,pEnd))
            return true;
    }
    return false;
}

/** \brief Verifies if a feature of tyhe specified type specifically overlaps the given range. Only valide for positions
 *
 * \param pStart long int: Begining of range
 * \param pEnd long int: End of range
 * \param pType featureType: Type of feature to check
 * \return bool: True if at least one overlap found
 *
 */
bool uGene::isOverlappingFeature(long int pStart, long int pEnd, featureType pType)
{
    for(auto & feature:m_featureVector)
    {
        if ( feature.getType()==pType && utility::isOverlap(feature.getStart(),feature.getEnd(),pStart,pEnd))
            return true;
    }
    return false;

}

/** \brief Utility function, returns the type of feature associated with a given string. OTHER returned if no match
 *
 * \param pType const std::string& : String to validate
 * \return featureType Type of feature mapped to it
 *
 */
featureType mapFeature(const std::string & pType )
{


    if (featureMap.count(  pType  ))
        return featureMap.find(pType)->second;
    else
        return featureType::OTHER;

}

/** \brief Utility function, return the string key associated with a feature type
 *
 * \param pType const featureType&: Feature type to check
 * \return std::string: associated string key
 *
 */
std::string featureString(const featureType& pType)
{
    for (auto & curPar:featureMap)
    {
        if (curPar.second==pType)
            return curPar.first;
    }
    return "";
}

unsigned long int uGene::featureCount(const featureType & pFeature)const
{
    if (pFeature==featureType::NULLFEATURE)
    {
        return m_featureVector.size();
    }
    else
    {
        return count_if(m_featureVector.begin(),m_featureVector.end(),[&](const uFeature & curfeature)
        {
            return (curfeature.getType()==pFeature);
        } );
    }
}

/**< uGeneChrom */

/**< Constructors */
uGeneChrom::uGeneChrom(const std::vector<uGene> & copyVec)
{
try {
    if (copyVec.size())
        setChr(copyVec.at(0).getChr());
    for (uGene elem: copyVec)
        addData(elem);
}
catch(...){throw;}
}


uGeneChrom::uGeneChrom(const std::string & pChrom, const std::vector<uGene> & pCopyVec)
{
    setChr(pChrom);
    for (uGene elem: pCopyVec)
        addData(elem);
}


/** \brief Total number of features, without counting genes
 *
 * \return unsigned long int
 *
 */
unsigned long int uGeneChrom::featureCount(const featureType & pFeature)const
{

    return accumulateSitesInfo([&](unsigned long int partialSum, const uGene & item) -> unsigned long int
    {
        return partialSum + item.featureCount(pFeature);
    }, 0ULL);
}



/** \brief Special function to add Data to a chrom.
 *
 *  As per uGene specified, an added item that shares an already used ID/Transcript combination will not be added as a new uGene but as a new feature
 *  of the previously existing item.
 *
 * \param newSite const uGene&: Item to add
 * \return void
 *
 */
void uGeneChrom::addData(const uGene& newSite)
try
{
    if (newSite.getChr()!=this->chr)
        throw ugene_exception_base()<<string_error("adding base to Chrom with non-matching scaffold/chr name");
    /**< Check for ID/Transcript combination if an ID is set*/
    if (newSite.getID().size())
    {
        /**< If we are trying to add but the existing ID/Transcript combination is set, fail */
        if (this->getIDCount(newSite.getID(),newSite.getTranscript()))
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
void uGeneChrom::addData(const uToken & pToken)
{
    /**< We first validate if there is a ID and potentially a Transcript associated. If not, treat it as normal */
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
            mainItr= this->_findGeneItr(tokID,tokTranscript);
        }
        for (; i<pToken.paramCount(token_param::START_POS); i++)
        {
            StrandDir dir=StrandDir::FORWARD;
            /**< Others as supplementary feature */
            if ( pToken.isParamSet(token_param::STRAND,i) && pToken.getParam(token_param::STRAND,i)=="-")
                dir=StrandDir::REVERSE;

            std::string featureID="", featureClass="", myOffset="0";

            if (pToken.isParamSet(token_param::FEATURE_ID,i))
                featureID=pToken.getParam(token_param::FEATURE_ID,i);

            if (pToken.isParamSet(token_param::PHASE,i))
                myOffset=pToken.getParam(token_param::PHASE,i);

            mainItr->addFeature(std::stoll(pToken.getParam(token_param::START_POS,i)),std::stoll(pToken.getParam(token_param::END_POS,i)),dir,mapFeature(pToken.getParam(token_param::FEATURE_TYPE,i)),featureID,std::stoi(myOffset));
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


/** \brief Get the number of genes that match the given ID. optionally, the number of genes that matched ID/Transcription combination (0 or 1 )
 *
 * \param pId const std::string&: ID to check
 * \param pTranscript const std::string&: Transcript to check
 * \return long int
 *
 */
long int uGeneChrom::getIDCount(const std::string & pId, const std::string & pTranscript)const
{

    if (pTranscript.size()==0){
            auto IdEqualFunct=[&](const uGene & item)
            {
                return (item.getID()==pId) ;
            };
            return  ( std::count_if(std::begin(VecSites), std::end(VecSites), IdEqualFunct) );
        }
    else
        {
        auto bothIdEqualFunct=[&](const uGene & item)
        {
           return ( item.getID()==pId && item.getTranscript()==pTranscript) ;
        };

        return  ( std::count_if(std::begin(VecSites), std::end(VecSites),bothIdEqualFunct ));
    }
}

/** \brief Return an iterator to the gene matching the ID/Transcript combination. return end() if none found. If "" passed, will return .end
 *
 * \param pId const std::string&: GeneID to match
 * \param pTranscript const std::string&: TranscriptID to match
 * \return typename std::vector<uGene>::iterator :Iterator to found element
 *
 */
typename std::vector<uGene>::const_iterator uGeneChrom::findGene(const std::string & pId, const std::string pTranscript)const
{
    if (pTranscript=="")
    {

        if (pId=="")
            return std::end(VecSites);

        return std::find_if(std::begin(VecSites),std::end(VecSites),[&](const uGene & item)
        {
            return(item.getID()==pId  ) ;
        });
    }
    else
    {

        return std::find_if(std::begin(VecSites),std::end(VecSites),[&](const uGene & item)
        {
            return(item.getID()==pId && item.getTranscript()== pTranscript ) ;
        });
    }
}

typename std::vector<uGene>::iterator uGeneChrom::_findGeneItr(const std::string & pId, const std::string pTranscript)
{
    if (pTranscript=="")
    {
        return std::find_if(std::begin(VecSites),std::end(VecSites),[&](uGene & item)
        {
            return(item.getID()==pId ) ;
        });
    }
    else
    {

        return std::find_if(std::begin(VecSites),std::end(VecSites),[&](uGene & item)
        {
            return(item.getID()==pId && item.getTranscript()== pTranscript ) ;
        });
    }
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


/** \brief Return a copy of the uGeneChrom structur
 *
 * \return uGeneChrom Copy of this
 *
 */
uGeneChrom uGeneChrom::getCopy()const
{
    uGeneChrom returnCopy = *this;
    return returnCopy;
}


/** \brief Similair to findNext, but will only find an item that contains a feature matching pType. Still works on current sort.
 *
 * \param pPosition long int  Position to from
 * \param pType featureType Feature type to filter on
 * \return typename std::vector<uGene>::const_iterator Iterator to the item, returns end() if not found.
 *
 */
typename std::vector<uGene>::const_iterator uGeneChrom::findNextWithFeature(long int pPosition, featureType pType)const
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
        if (upper==VecSites.end()){
            return VecSites.end();
            }

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
        std::cerr << "Calling findNextGeneWithFeature and you did not provide an appropriate get function" <<std::endl;
#endif
        throw e;
    }
    catch (...)
    {
#ifdef DEBUG
        std::cerr << "Catch and re-throw from findNextGeneWithFeature" <<std::endl;
#endif
        throw ;
    }
}


/** \brief Similair to findPreceding, but will only find an item that contains a feature matching pType. Still works on current sort.
 *
 * \param pPosition long int  Position to from
 * \param pType featureType Feature type to filter on
 * \return typename std::vector<uGene>::const_iterator Iterator to the item, returns end() if not found.
 *
 */
//typename std::vector<uGene>::const_iterator uGeneChrom::findPrecedingWithFeature(long int pPosition, featureType pType)const
//{
//    try
//    {
//        /**< If unsorted, fail */
//        if (VecSites.size()==0)
//            return VecSites.end();
//        if (m_isSorted==false)
//            throw unsorted_throw() <<string_error("findPrecedingSite called on unsorted vector \n") ;
//        if ((sortGetStart==nullptr)||(sortGetEnd==nullptr))
//            throw ugene_exception_base() <<string_error(" findPrecedingSite called on chrom without appropriate start or end function\n") ;
//        auto comp = [&] (const uGene &item1, const float &pPos)
//        {
//            return ( sortGetStart(&item1)< pPos && item1.hasFeatureType(pType) );
//        };
//
//        /**< Compare, sort Value */
//        auto lower = std::lower_bound(VecSites.begin(), VecSites.end(), pPosition, comp);
//        /**< If result is our first item, then no item precedes it */
//        if ((lower==VecSites.begin())){
//            return VecSites.end();
//        }
//        /**< If result is end, every item precedes the value */
//        /**<Return item precedes and as such is LESS then position. if no item was found, last item is closest to value  */
//        lower--;
//        return (lower);
//    }
//    catch (unsorted_throw & e)
//    {
//#ifdef DEBUG
//        std::cerr << "FindPrecedingSite called on unsorted vector" <<std::endl;
//#endif
//        throw e;
//    }
//    catch (ugene_exception_base & e)
//    {
//#ifdef DEBUG
//        std::cerr << "Calling findPrecedingSite and you did not provide an aproriate get function" <<std::endl;
//#endif
//        throw e;
//    }
//}



//TODO Complete this
/** \brief Greedy search as our data is not sorted. In fact, we need to pass every element to check if it's boundary overlaps. Will search passed on current sort
 *
 * \param pPosition long int  Position to from
 * \param pType featureType Feature type to filter on
 * \return typename std::vector<uGene>::const_iterator Iterator to the item, returns end() if not found.
 *
 */
//std::pair<typename std::vector<uGene>::const_iterator,typename std::vector<uFeature>::const_iterator> uGeneChrom::findNextFeature(long int pPosition, featureType pType)const
//{
//    typename std::vector<uGene>::const_iterator curGene=VecSites.end();
//    typename std::vector<uFeature>::const_iterator curFeature;
//
//    long int curPos=0;
//    for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
//    {
//        if (it->hasFeatureType(pType) && it->getBoundaryEnd() > pPosition)
//        {
//            for (auto featureItr= it->featureBegin(); featureItr!=it->featureEnd(); featureItr++)
//            {
//                if (featureItr->getType()==pType && featureItr->getStart()>pPosition && featureItr->getStart()>curPos)
//                {
//                    curGene=it;
//                    curPos= featureItr->getStart();
//                }
//
//            }
//        }
//    }
//    return (std::make_pair(curGene,curFeature));
//}



/**< uGeneExperiment */
/** \brief Wrapper around the Chrom level function
 *
 * \param pChr std::string Scaffold to run the function on
 * \param pPosition long int Position to search from
 * \param pType featureType FeatureType to find
 * \return typename std::vector<uGene>::const_iterator Iterator to gene
 * \sa findNextGeneWithFeature
 */
typename std::vector<uGene>::const_iterator uGeneExperiment::findNextWithFeature(std::string pChr, long int pPosition, featureType pType)const
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
/** \brief Wrapper around the Chrom level function
 *
 * \param pChr std::string Scaffold to run the function on
 * \param pPosition long int Position to search from
 * \param pType featureType FeatureType to find
 * \return typename std::vector<uGene>::const_iterator Iterator to gene
 * \sa findPrecedingGeneWithFeature
 */
//typename std::vector<uGene>::const_iterator uGeneExperiment::findPrecedingGeneWithFeature(std::string pChr,long int pPosition, featureType pType)const
//{
//    try
//    {
//        if (!ExpMap.count(pChr))
//        {
//            throw param_throw()<<string_error("Failling in uGenericNGSExperiment::findNextSite, value "+pChr+" does not exist.\n");
//        }
//        auto tempChrom = getpChrom(pChr);
//        return tempChrom->findPrecedingWithFeature(pPosition,pType);
//    }
//    catch(...)
//    {
//        throw;
//    }
//}

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

/** \brief Return number of feature in the exp, ignoring the main feature. Optionally, only count certain type of features
 *
 * \param pFeature=featureType::NULLFEATURE const featureType&
 * \return unsigned int
 *
 */
unsigned long int uGeneExperiment::featureCount(const featureType &pFeature)const
{
    return accumulateChromsInfo([&](unsigned long int partialSum, const uGeneChrom & item) -> unsigned long int
    {
        return partialSum + item.featureCount(pFeature);
    }, 0ULL);
}

/** \brief Wrapper around the chrom level function, adding a Token to uGene is mildly more complicated.
 *
 * \param pToken uToken& uToken to parse and add
 * \return void
 *
 */
void uGeneExperiment::addData(const uToken & pToken)
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
