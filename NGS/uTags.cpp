#include <stdio.h>
#include <string.h>
#include <sstream>
#include "uTags.h"
#include "uBasicNGS.h"
#include "uRegion.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
namespace NGS
{

using namespace std;
/** \brief Default constructor, not PE and positive strand
 */
uTags::uTags():uGenericNGS()
{}

/**< From uTokens */
uTags::uTags(uToken pToken)try:
    uGenericNGS(pToken)
{
    if (pToken.isParamSet(token_param::SEQ_NAME))
        setName(pToken.getParam(token_param::SEQ_NAME));
    if (pToken.isParamSet(token_param::SEQUENCE))
        setSequence(pToken.getParam(token_param::SEQUENCE));
    if (pToken.isParamSet(token_param::CIGAR))
        setCigar(pToken.getParam(token_param::CIGAR));
    if (pToken.isParamSet(token_param::MAP_SCORE))
        setMapQual(utility::stoi(pToken.getParam(token_param::MAP_SCORE)));
    if (pToken.isParamSet(token_param::PHRED_SCORE))
        setPhred(pToken.getParam(token_param::PHRED_SCORE));
    if (pToken.isParamSet(token_param::FLAGS))
        setFlag(utility::stoi(pToken.getParam(token_param::FLAGS)));
    if (pToken.isParamSet(token_param::TEMPLATE_LENGHT))
        setPELenght(utility::stoi(pToken.getParam(token_param::TEMPLATE_LENGHT)));

}
catch(ugene_exception_base &e)
{
#ifdef DEBUG
    std::cerr << "Error in uTags(uToken)." <<std::endl;
#endif
    e<<tag_error(*this);
    throw e;
}


/** \brief Copy constructor, with init list
 * \param otherItem: uBasicNGS  object
 */
uTags::uTags(const uBasicNGS & otherItem):uGenericNGS(otherItem.getChr(),otherItem.getStart(),otherItem.getEnd(),otherItem.getStrand()),name(nullptr),phredScore(nullptr),m_cigar(nullptr)
{
    this->setScoreVector(otherItem.getScoreVector());
}

/** \brief Copy constructor, with init list
 * \param otherItem: uRegion  object
 */
uTags::uTags(const uRegion & otherItem):uGenericNGS(otherItem.getChr(),otherItem.getStart(),otherItem.getEnd(),otherItem.getStrand()),name(nullptr),phredScore(nullptr),m_cigar(nullptr)
{
    this->setScoreVector(otherItem.getScoreVector());
}

/** \brief Default constructor with init list, implicitly sets strand
 *
 * \param chr: name of the chromosome
 * \param start: beginning position of the tag
 * \param end: ending position of the tag
 */
uTags::uTags(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand):name(nullptr),phredScore(nullptr),m_cigar(nullptr)
{
    try
    {
        setStartEnd(pStart,pEnd);
        setChr(pChr);
        setStrand(pStrand);
    }
    catch(elem_throw & e)
    {
#ifdef DEBUG
        cerr << "Throwing in uTags constructor" <<endl;
#endif

        string trace;
        if (std::string const * ste =boost::get_error_info<string_error>(e) )
            trace=*ste;

        e << string_error(trace+"Failling in uTags constructor, parameters are"+pChr+" "+std::to_string(pStart)+" "+std::to_string(pEnd)+"\n");

        throw e;
    }
}

uTags::uTags(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore)try :
    uGenericNGS(pChr,pStart,pEnd,pStrand,pScore),name(nullptr),phredScore(nullptr),m_cigar(nullptr)
{
    /* try
     {
         setStartEnd(pStart,pEnd);
         setChr(pChr);
         setStrand(pStrand);
         setScore(pScore);
     }*/
}
catch(construct_elem_throw & e)
{
#ifdef DEBUG
    cerr << "Throwing in uTags constructor" <<endl;
#endif
    string trace;
    addStringError(e,("Failling in uTags constructor, parameters are"+pChr+" "+std::to_string(pStart)+" "+std::to_string(pEnd)+"\n") );
    e << tag_error(*this);
    throw e;
}

uTags::uTags(std::string pChr, long long int pStart, long long int pEnd, float pScore)try :
    uGenericNGS(pChr,pStart,pEnd,pScore),name(nullptr),phredScore(nullptr),m_cigar(nullptr)
{}
catch(construct_elem_throw &e)
{
#ifdef DEBUG
    cerr << "Throwing in uTags constructor" <<endl;
#endif
    addStringError(e,"Throwing in uTags(string,long long int, long long int,float)");
    e << tag_error(*this);
    throw e;
}


/** \brief Default destructor
 */
uTags::~uTags()
{

    if (name!=nullptr)
    {
        delete []name;
        name = nullptr;
    }

    if (phredScore!=nullptr)
    {
        delete []phredScore;
        phredScore = nullptr;
    }

    if (m_cigar!=nullptr)
    {
        delete []m_cigar;
        m_cigar = nullptr;
    }
}

//TODO TEST THIS
/** \brief Copy Constructor for UTags, necessary due to our char*
 * \param UTags& const copy_from: tag to copy
 */
uTags::uTags(const uTags& copy_from):uGenericNGS(copy_from),name(nullptr),phredScore(nullptr),m_cigar(nullptr)
{
    try
    {
        if (copy_from.name!=nullptr)
        {
            name = new char[(strlen(copy_from.name)+1)];
            strcpy(name,copy_from.name);
        }

        if (copy_from.m_cigar!=nullptr)
        {
            m_cigar = new char[(strlen(copy_from.m_cigar)+1)];
            strcpy(m_cigar,copy_from.m_cigar);
        }

        if (copy_from.phredScore!=nullptr)
        {
            phredScore = new char[(strlen(copy_from.phredScore)+1)];
            strcpy(phredScore,copy_from.phredScore);
        }
        m_strand=copy_from.m_strand;
        PELenght=copy_from.PELenght;
        sequence=copy_from.sequence;
        mapScore=copy_from.mapScore;
        Unmapped= copy_from.Unmapped;
    }
    catch(elem_throw & e)
    {
#ifdef DEBUG
        cerr << " Failling in copy constructor, uTags(uTags&) " <<endl;
#endif
        string trace;
        if (std::string const * ste =boost::get_error_info<string_error>(e) )
            trace=*ste;
        e << string_error(trace+"Failling in copy constructor, uTags(uTags&) \n");
        e<< tag_error(copy_from);
        throw e;
    }
    catch(std::exception &e)
    {
        throw e;
    };
}

/** \brief Overloaded assignement operator
 * \param uTags const& source
 * \return uTags& destination
 */
uTags& uTags::operator= (uTags const& assign_from)
{
    if (this == &assign_from) return *this;

    uGenericNGS::operator= (assign_from);

    if (name!=nullptr)
    {
        delete []name;
        name = nullptr;
    }

    if (phredScore!=nullptr)
    {
        delete []phredScore;
        phredScore = nullptr;
    }

    if (m_cigar!=nullptr)
    {
        delete []m_cigar;
        m_cigar = nullptr;
    }

    if (assign_from.name!=nullptr)
    {
        name = new char[(strlen(assign_from.name)+1)];
        strcpy(name,assign_from.name);
    }

    if (assign_from.m_cigar!=nullptr)
    {
        m_cigar = new char[(strlen(assign_from.m_cigar)+1)];
        strcpy(m_cigar,assign_from.m_cigar);
    }

    if (assign_from.phredScore!=nullptr)
    {
        phredScore = new char[(strlen(assign_from.phredScore)+1)];
        strcpy(phredScore,assign_from.phredScore);
    }
    m_strand=assign_from.m_strand;
    PELenght=assign_from.PELenght;
    sequence=assign_from.sequence;
    mapScore=assign_from.mapScore;
    Unmapped= assign_from.Unmapped;
    return *this;
}


uTags uTags::getCopy() const
{
    uTags copyObj=*this;
    return copyObj;
}

bool uTags::isEqual(const uTags & pCompared) const
{
    return ((this->getChr()==pCompared.getChr())&&
            (this->getStrand()==pCompared.getStrand())&&
            (this->getStart()==pCompared.getStart())&&
            (this->getEnd()==pCompared.getEnd())&&
            (this->getScoreVector()==pCompared.getScoreVector())&&
            (this->getMapQual()==pCompared.getMapQual())&&
            (this->getSequence()==pCompared.getSequence())&&
            (this->getName()==pCompared.getName())&&
            (this->getPhred()==pCompared.getPhred())&&
            (this->getCigar()==pCompared.getCigar())&&
            (this->isMapped()==pCompared.isMapped())&&
            (this->getFlag()==pCompared.getFlag())&&
            (this->getPeLenght()==pCompared.getPeLenght()));

}

/** \brief Returns true if the tag is unmapped or incorrectly mapped
 *
 * \return bool
 *
 */
bool uTags::isMapped() const
{
    return Unmapped;
}
/** \brief Set the mapped value of the tag.
 *
 * \param pmapped bool : Value to set the status to
 * \return void
 *
 */
void uTags::setMapped(bool pmapped)
{
    Unmapped=pmapped;
}


/** \brief Set the cigar value
 *
 * \param pcigar std::string : String to set cigar to
 * \return void
 *
 */
void uTags::setCigar(std::string pCigar)
{

    //   if (!( utility::validateCigar(pCigar))
    //         throw param_throw()<<string_error("Invalide m_cigar flag passed to SetCigar. Flag is: "+pCigar+"\n" )
    if (m_cigar!=nullptr)
    {
        delete []m_cigar;
        m_cigar = nullptr;
    }
    try
    {
        m_cigar = new char[pCigar.size()+1];
        int lenght=pCigar.copy(m_cigar,pCigar.size(),0 );
        m_cigar[lenght]='\0';

    }
    catch(std::exception & e)
    {
#ifdef DEBUG
        std::cerr <<"Failled allocation in setCigar()";
#endif
        throw  e;
    }
}

/** \brief get the Cigar string associated with the sequence
 *
 * \return std::string
 *
 */
std::string uTags::getCigar() const
{
    std::string returnStr;
    if (m_cigar==nullptr)
        returnStr="";
    else
        returnStr=m_cigar;
    return returnStr;
}

/** \brief Set the Sam Flag associated with the sequence
 *
 * \param pflag int : The sam flag to set
 * \return void
 *
 */
void uTags::setFlag(int pflag)
{
    flag=pflag;
}
/** \brief Get the Sam flag associated with the sequence.
 *
 * \return int : The Sam flag associated.
 *
 */
int uTags::getFlag() const
{
    return flag;
}

/** \brief Set the sequence associated with the element
 *
 *  Sets the sequence associated with the element. The sequence needs to be either
 *  null ("") or of size equal to the element.
 * \exception param_throw : Thrown when parameter size neighter null or equal to element getLenght()
 * \param pSeq std::string : The sequence to set
 * \return void
 *
 */
void uTags::setSequence(std::string pSeq)
{
    if(((int)pSeq.size()!=0)&&((int)pSeq.size()!=getLenght()))
        throw param_throw()<< string_error("Failling in seqSequence. Sequence size neither null or equal to element size.");
    sequence=pSeq;
}
/** \brief Get the sequence associated with the element.
 *
 * \return std::string : The sequence associated with the element.
 *
 */
std::string uTags::getSequence() const
{
    return sequence;
}

/** \brief Set the PhredScore assocaited with the sequence
 *
 * \param Phred std::string : The sequence to set
 * \exception
 * \return void
 *
 */
void uTags::setPhred(std::string Phred)
{
    if (phredScore!=nullptr)
    {
        delete []phredScore;
        phredScore = nullptr;
    }

    try
    {
        phredScore = new char[Phred.size()+1];
        int lenght=Phred.copy(phredScore,Phred.size(),0 );
        phredScore[lenght]='\0';
    }
    catch(std::exception & e)
    {
#ifdef DEBUG
        std::cerr <<"Failled allocation in setPhred()";
#endif
        throw e;
    }


}

/** \brief get the PhredScore associated with the sequence
 *
 * \return std::string : The PhredScore vector
 *
 */
std::string uTags::getPhred() const
{
    std::string returnStr;
    if (phredScore==nullptr)
        returnStr=="";
    else
        returnStr=phredScore;
    return returnStr;
}

/** \brief Set the ID associated with the element.
 *
 *  This sets the ID associated with the element. There are no conditions attached to the ID
 *
 * \param pName std::string : ID to set the name to
 * \exception Thrown when new allocations fails.
 * \return void
 *
 */
void uTags::setName(std::string pName)
{
    if (name!=nullptr)
    {
        delete []name;
        name = nullptr;
    }

    try
    {
        name= new char [pName.size()+1];
        int lenght = pName.copy(name, pName.size(),0);
        name[lenght]='\0';
    }
    catch(std::exception & e)
    {
#ifdef DEBUG
        std::cerr <<"Failled allocation in setName()";
#endif

        throw e ;
    }

}

/** \brief Get the name/ID associated with the element
 *
 * \return std::string : The returned ID
 *
 */
std::string uTags::getName() const
{
    std::string returnStr;
    if (name==nullptr)
        returnStr="";
    else
        returnStr=name;
    return returnStr;
}

/** \brief True if the paired end lenght is above 0
 *
 * \return bool True if PElenght is set
 *
 */
bool uTags::isPE() const
{
    return PELenght;
}

/** \brief
 *
 * \param lenght int : Value to set PELenght to.
 * \exception : param_throw(): Throw if parameter is < 0.
 * \return void
 *
 */
void uTags::setPELenght(int lenght)
{

    if  (lenght <0)
        throw param_throw()<<string_error("Throwing in setPELenght. Set an invalid PE lenght<0");
    PELenght=lenght;

}
/** \brief Return the PELenght of the element
 *
 * \return int
 *
 */
int uTags::getPeLenght() const
{
    return PELenght;
}

/** \brief Set the mapping quality of the element, max value is 255
 *
 * \param score short int : Value to set MapQuality to
 * \return void
 *
 */
void uTags::setMapQual(short int score)
{
    mapScore=score;
}

/** \brief Return the mapping quality of the element.
 *
 * \return short int: Mapping value
 *
 */
short int uTags::getMapQual() const
{
    return mapScore;
}

/** \brief Returns an element that represents the "completed" tag of a paired end item.
 *
 * \return uTags The completed tag. If no PE data, returns *this
 *
 */
uTags uTags::getCompletedCopy()const
{
    try
    {
        if (PELenght==0)
            return this->getCopy();
        else
        {
            uTags rtnCopy= this->getCopy();
            if (rtnCopy.getStrand()==StrandDir::FORWARD)
                rtnCopy.setEnd(rtnCopy.getEnd()+rtnCopy.getPeLenght());
            else
                rtnCopy.setStart(rtnCopy.getStart()-rtnCopy.getPeLenght());
            return rtnCopy;
        }
    }
    catch(ugene_exception_base & e)
    {
        addStringError(e, "Failling in getCompletedTag, completing tag would invalid start or end position");
        e <<tag_error(*this);
        throw e;
    }
}

/** \brief Prints a human readable version of the element in no particular format.
 *
 * \param pOut std::ostream& Output to write to.
 * \return void
 *
 */
void uTags::print(std::ostream &pOut) const
{
    pOut<<"Chrom: "<<getChr()<<std::endl;
    pOut<<"Start: "<<utility::to_string(getStart())<<std::endl;
    pOut<<"End: " <<utility::to_string(getEnd())<<std::endl;
    if (Unmapped)
    {
        pOut<<"Tag is not mapped to reference or mapped with errors"<<std::endl;
    }
    if (m_score.size()>0)
    {
        pOut<<"Scores: ";
        for(auto value: m_score )
            pOut<<value<<" ";
        pOut<<std::endl;
    }
    if (name!=nullptr)
    {
        pOut<<"Name: " <<name<<std::endl;
    }
    if (phredScore!=nullptr)
    {
        pOut<<"Pred Score: " <<phredScore<<std::endl;
    }

    if (m_cigar!=nullptr)
    {
        pOut<<"cigar: " <<m_cigar<<std::endl;
    }
    if (this->isPE())
    {
        pOut<<"Is paired, PE lenght is " <<utility::to_string(getPeLenght())<<std::endl;
    }
    else
    {
        pOut<<"Is not paired "<<std::endl;
    }
    pOut<<"Map score: " <<utility::to_string(getMapQual())<<std::endl;
    pOut<<"Flag: " <<utility::to_string(getFlag())<<std::endl;
    pOut<<"Seq: " <<getSequence()<<std::endl;

}

/** \brief Create the parser Token associated with the uTag
 *
 * \return uToken The token returned
 *
 */
uToken uTags::createToken() const
{
    std::stringstream ss;
    ss << "CHR\t"<<this->getChr()<<"\nSTART_POS\t"<<this->getStart()<<"\n" << "END_POS\t"<<this->getEnd()<<"\n";
    if (getScoreCount()>0)
    {
        ss << "SCORE\t"<<this->getScore()<<"\n";
    }
    if (getStrand()==StrandDir::FORWARD)
        ss << "STRAND\t"<<"+"<<"\n";
    else
        ss << "STRAND\t"<<"-"<<"\n";

    ss << "MAP_SCORE\t"<<utility::to_string(this->getMapQual())<<"\nFLAGS\t"<<utility::to_string(this->getFlag())<<"\n";
    if (getSequence()!="")
        ss << "SEQUENCE\t"<<getSequence()<<"\n";

    if (name!=nullptr)
    {
        ss<<"SEQ_NAME\t"<<name<<"\n";
    }
    if (phredScore!=nullptr)
    {
        ss<<"PHRED_SCORE\t"<<phredScore<<"\n";
    }
    if (m_cigar!=nullptr)
    {
        ss<<"CIGAR\t" <<m_cigar<<"\n";
    }
    try
    {
        return uToken(ss);
    }
    catch(uToken_exception_base &e)
    {
        addStringError(e, "Failed while creating token in uGenericNGS::createToken()");
        throw e;
    }

}


/** \brief Copy constructor
 *
 * \param copyCop : The object to instaciate  from.
 */
uTagsChrom::uTagsChrom(const uGenericNGSChrom<uTagsChrom,uTags> & copyCop)
{
    VecSites=copyCop.returnVecData();
    chr= copyCop.getChr();
    m_isSorted=copyCop.getSortedStatus();
    sortGetStart=copyCop.getStartFunct();
    sortGetEnd=copyCop.getEndFunct();
    m_comptFunc=copyCop.getCompFunct();
    chromSize=copyCop.getChromSize();
}

uTagsChrom::uTagsChrom(const uTagsChrom& initFrom)
{

    VecSites=initFrom.returnVecData();
    chr= initFrom.getChr();
    m_isSorted=initFrom.m_isSorted;
    sortGetStart=initFrom.sortGetStart;
    sortGetEnd=initFrom.sortGetEnd;
    m_comptFunc=initFrom.m_comptFunc;
    chromSize=initFrom.chromSize;

}

uTagsChrom& uTagsChrom::operator=(const uTagsChrom& copFrom)
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

uTagsChrom uTagsChrom::getCopy() const
{
    uTagsChrom copyObj=*this;
    return copyObj;
}

uTagsChrom::uTagsChrom(const uRegionChrom & pCopyChrom)
{
    setChr(pCopyChrom.getChr());
    chromSize=pCopyChrom.getChromSize();
    for (auto itr= pCopyChrom.begin(); itr!=pCopyChrom.end(); itr++  )
        addData(uTags(*itr));

}

uTagsChrom::uTagsChrom(const uBasicNGSChrom & pCopyChrom)
{
    setChr(pCopyChrom.getChr());
    chromSize=pCopyChrom.getChromSize();
    for (auto itr= pCopyChrom.begin(); itr!=pCopyChrom.end(); itr++  )
        addData(uTags(*itr));
}

//TODO Re-write this, same comment as Experiment level function
/** \brief For a given start and end on this Chromosome, return the it's continuous density signal
 * \doto Use our standard overlap function
 *
 * \param start : Begin of our region
 * \param end : End of our regions
 * \param overlap : boolean, do we expect overlap or englobed
 *
 * \return
 *
 */
std::vector<float> uTagsChrom::getRegionSignal(int start, int end, bool overlap)
{
    vector<float> tempSignal;
    int signalStart, signalEnd, signalSize;

    tempSignal.resize((end-start)+1);
    std::vector<uTags>::iterator iterVec;


    /**< Assuming the data is sorted, this could be heavily optimized. */
    int pos=0;
    iterVec=VecSites.begin();

    for (iterVec=(iterVec+pos) ; iterVec != VecSites.end(); ++iterVec)
    {
        /**< If our start is passed our region  we stop. time saver. */
        if (iterVec->getStart()> end)
            break;

        if (overlap==false)
        {
            if (utility::isRegionAInsideRegionB(iterVec->getStart(), iterVec->getEnd(), start, end))
            {
                /**< Increase signal for every overlapping position of the tag. */
                signalStart= (iterVec->getStart()-start);
                signalEnd= (iterVec->getEnd()-start);
                signalSize=(signalEnd-signalStart  );
                for (int i=0; i<=signalSize; i++ )
                {
                    tempSignal.at(signalStart+i)++;
                }
            }
        }
        else if (utility::checkOverlap(iterVec->getStart(), iterVec->getEnd(), start, end))
        {
            signalStart= (iterVec->getStart()-start);
            signalEnd= (iterVec->getEnd()-start);
            if (signalStart<0)
                signalStart=0;

            if (signalEnd>(end-start))
                signalEnd= (end-start);

            signalSize=(signalEnd-signalStart  );

            for (int i=0; i<=signalSize; i++ )
            {
                tempSignal.at(signalStart+i)++;
            }

        }
    }
    return tempSignal;
}

template <class _OTHER_>
uTags uTagsChrom::generateRandomSite(const int size_,std::mt19937& engine,const _OTHER_ &exclList, const int sigma, const std::string ID) const
{
    //TODO Sanity check here to make sure it is possible to generate the asked for tag.
    uTags returnTag;

    bool found=false;
    int size = size_;

    int max = this->getChromSize();

    while (!found)
    {
        {
            uTags temptag;
            if (sigma!=0)
            {
                std::normal_distribution<float> gaussian(size, sigma);
                size = (int)gaussian(engine);
            }

            if (size>=1)
            {
                int shift = size/2;
                //Generating our distribution at each call is probably needlesly heavy.. check to optimize this in time.
                std::uniform_int_distribution<int> unif((shift+1), (max-shift));
                int center = unif(engine);
                temptag.setEnd(center+shift);
                temptag.setStart(center-shift);
                if ((exclList.getSubset(temptag.getStart(),temptag.getEnd())).count()==0)
                {
                    found=true;
                    returnTag=temptag;
                    returnTag.setChr(this->getChr());
                    returnTag.setName(ID);
                }
            }
        }
    }

    return returnTag;
}

//TODO Re-write this to be a templated free function. While densignal signal makes more sense for Tags, they could be used from any structure
/** \brief Generate our density signal for a given region and chrom
 *
 * \param chrom std::string : The specified chrom
 * \param start int : start position of the region
 * \param end int : end position of the region
 * \param overlap bool : Type of overlap
 * \return std::vector<float>
 *
 */
std::vector<float> uTagsExperiment::getRegionSignal(std::string chrom, int start, int end, bool overlap)
{
    vector<float> returnVec;

    std::map<std::string,uTagsChrom>::iterator iterMap;
    uTagsChrom* ptempChrom;

    ptempChrom=&(ExpMap[chrom]);

    returnVec= ptempChrom->getRegionSignal(start, end, overlap);

    return returnVec;
}


/** \brief Return a copy of the experiment
 *
 * \return uTagsExperiment
 *
 */
uTagsExperiment uTagsExperiment::getCopy() const
{
    uTagsExperiment copyObj=*this;
    return copyObj;
}



/** \brief Load the next blockSize number of entries from the given Bamreader. if not more data available, stop reader
 *
 * \param pReader BamReader& BamReader to use
 * \param blockSize=1; int : Number of entries to read
 * \param pLoadCore bool : If true, skips certain string data ( sequence, phred scores, Name )
 * \return uTagsExperiment
 *
 */
void uTagsExperiment::loadWithBamTools(BamTools::BamReader& pReader, int pBlockSize, bool pLoadCore )
{
    try {
        if (pLoadCore)
            loadWithBamTools_Core(pReader, pBlockSize);
        else
            loadWithBamTools_All(pReader, pBlockSize);
    }catch(...){
    throw;}

}

void uTagsExperiment::loadWithBamTools_Core(BamTools::BamReader& pReader, int pBlockSize)
{
 try {
   int countLoaded=0;
    BamTools::BamAlignment m_BufferAlignement;
    /**< if no buffer, load data */
    while(countLoaded<pBlockSize)
    {
        if (pReader.GetNextAlignmentCore(m_BufferAlignement))
        {
            uTags tagToAdd(pReader.GetReferenceData().at(m_BufferAlignement.RefID).RefName,(m_BufferAlignement.Position+1),m_BufferAlignement.GetEndPosition());

            if (m_BufferAlignement.IsReverseStrand())
                tagToAdd.setStrand(NGS::StrandDir::REVERSE);

            tagToAdd.setFlag(m_BufferAlignement.AlignmentFlag);
            tagToAdd.setMapQual(m_BufferAlignement.MapQuality);
            tagToAdd.setPELenght(m_BufferAlignement.InsertSize);

            std::string cigar;
            for(BamTools::CigarOp & cigarItem: m_BufferAlignement.CigarData)
            {
                cigar+= ( std::to_string(cigarItem.Type)+std::to_string(cigarItem.Length));
            }
            tagToAdd.setCigar(cigar);
            this->addData(tagToAdd);
            countLoaded++;
        }
        else
            break;
    }
}
catch(...){

throw;
}
}
void uTagsExperiment::loadWithBamTools_All(BamTools::BamReader& pReader, int pBlockSize)
{
    try {
  int countLoaded=0;
    BamTools::BamAlignment m_BufferAlignement;
    /**< if no buffer, load data */
    while(countLoaded<pBlockSize)
    {
        if (pReader.GetNextAlignment(m_BufferAlignement))
        {
            /**< Bam is 0 based, SAM 1 based. To map between them, +1 to Bam start positions. */
            uTags tagToAdd(pReader.GetReferenceData().at(m_BufferAlignement.RefID).RefName,(m_BufferAlignement.Position+1),(m_BufferAlignement.GetEndPosition()));

            if (m_BufferAlignement.IsReverseStrand())
                tagToAdd.setStrand(NGS::StrandDir::REVERSE);

            tagToAdd.setSequence(m_BufferAlignement.QueryBases);
            tagToAdd.setFlag(m_BufferAlignement.AlignmentFlag);
            tagToAdd.setMapQual(m_BufferAlignement.MapQuality);
            tagToAdd.setPELenght(m_BufferAlignement.InsertSize);
            tagToAdd.setPhred(m_BufferAlignement.Qualities);
            tagToAdd.setName(m_BufferAlignement.Name);

            std::string cigar;
            for(BamTools::CigarOp & cigarItem: m_BufferAlignement.CigarData)
            {
                cigar+= ( std::to_string(cigarItem.Type)+std::to_string(cigarItem.Length));
            }
            tagToAdd.setCigar(cigar);
            this->addData(tagToAdd);
            countLoaded++;
        }
        else{
            break;
        }
    }

}
catch(...)
{
    throw;
}

}
} // End of namespace NGS



/**< Note, this is deprecated and no longer used anywhere in the library.
    It is preserved "in-case" as in some very rare situations this could be of use.
 */
namespace factory
{
NGS::uTags makeTagfromSamString(std::string samString, bool minimal)
{
    std::stringstream Infostream;
    int size,ourFlag;
    Infostream.str(samString);
    std::string ourChr, ourStart, ourEnd, ourName;
    std::string phre, cig, seq;
    std::string temp,tempname;
    /**< Read name */
    try
    {
        if (!minimal)
        {
            Infostream>> tempname;
            if (tempname.find("/")!=std::string::npos)
            {
                tempname.erase(tempname.find("/"));
            }
            ourName=tempname;
        }
        else
            Infostream>>temp;

        /**< Read flag */
        Infostream>>ourFlag;
        Infostream>>ourChr;

        //returnTag.setChr(ourChr);

        Infostream>>ourStart;

        int mapScore;
        Infostream>>mapScore;
        /**< Cigar */
        Infostream>>cig;

        /**<  Parse Cigar here */
        /**< Find first letter , read number before, repeat */
        /**< Easier with Regex,  but what the hell */
        int curPos=0;
        std::string substr;
        size=0;
        for (unsigned int i=0; i< cig.size(); i++)
        {
            /**< If isAlpha then check previous numbers */
            if (isalpha(cig.at(i)))
            {
                /**< If a count value */
                char temp;
                temp = cig.at(i);
                if ((temp=='M')||(temp=='I')||(temp=='S')||(temp=='X')||(temp=='+'))
                {
                    substr= cig.substr(curPos, (i-curPos));
                    size+= atoi(substr.c_str());
                }
                curPos=(i+1);
            }
        }

//   ourEnd=(utility::stoi(ourStart)+(size-1));
        NGS::uTags returnTag(ourChr,utility::stoi(ourStart),utility::stoi(ourStart)+(size-1) );
        returnTag.setName(ourName);
        returnTag.setFlag(ourFlag);
        returnTag.setMapQual(mapScore);
        /**< name of next mate */
        Infostream>>temp;
        /**< Pos of next mate */
        Infostream>>temp;
        /**< Template lenght for PE readas */
        int PELenght;
        Infostream>>PELenght;
        PELenght= abs(PELenght);
        // cerr << "PeLenght is" <<
        /**< Sequence */
        Infostream>>seq;
        /**< Pred score of every position. */
        Infostream>>phre;
        /**< Sam flag */
        /**< StrongType this */
        if (ourFlag&0x10)
            returnTag.setStrand(NGS::StrandDir::REVERSE);// ='-';
        else
            returnTag.setStrand(NGS::StrandDir::FORWARD);// ='+'; */
        if (ourFlag&0x4)
            returnTag.setMapped(false);
        else
            returnTag.setMapped(true);

        /**< PE validation */
        /**< If PE and mate is aligned */
        if ((ourFlag&0x1)&&(ourFlag&0x2))
        {
            returnTag.setPELenght(PELenght);
        }
        else
        {
            returnTag.setPELenght(0);
        }
//    pMate = NULL;
        /**< if we want to keep, we store, otherwise scope will erase our data */
        if (!minimal)
        {
            returnTag.setPhred(phre) ;

            returnTag.setCigar(cig);
        }
        /**< Move semantics */
        return returnTag;

    }
    catch(NGS::elem_throw & e)
    {
        std::string trace;
#ifdef DEBUG
        cerr << "catching in factory on elem_throw" <<endl;
#endif

        if (std::string const * ste =boost::get_error_info<NGS::string_error>(e) )
            trace=*ste;

        e << NGS::string_error(trace+"Failling in factory:makeTagfromSamString constructor,  on string \n"+samString);
        throw e;
    }
    catch(...)
    {
#ifdef DEBUG
        cerr << "catching in factory on general throw" <<endl;
#endif

        NGS::elem_throw e;
        e << NGS::string_error("we threw in makeTagfromSamString trying the next string \n"+samString);
        throw e;
    }

}
}
