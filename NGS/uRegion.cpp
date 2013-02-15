#include "uFormats.h"
#include "uRegion.h"
#include "uTags.h"
#include "uBasicNGS.h"
#include "IO/uToken.h"

using namespace std;
namespace NGS {
/** \brief Default constructor
 */
uRegion::uRegion()
{}

/** \brief Region constructor, with start, end and chr
 *
 * \param ourchr std::string  : Chrom to set
 * \param ourstart int        : Start position
 * \param ourend int          : End position
 *
 */
uRegion::uRegion(std::string pChr, long long int pStart, long long int pEnd,StrandDir pStrand) try : uGenericNGS(pChr, pStart,pEnd,pStrand)
{}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(string,int int)");
    e << region_error(*this);
    throw e;
}
uRegion::uRegion(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore)try : uGenericNGS(pChr,pStart, pEnd,pStrand,pScore)
{}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(string,int int,Strandir,float)");
    e << region_error(*this);
    throw e;
}
uRegion::uRegion(std::string pChr, long long int pStart, long long int pEnd, float pScore )try : uGenericNGS(pChr,pStart,pEnd,pScore)
{}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(string,int, int,float)");
    e << region_error(*this);
    throw e;
}

uRegion::uRegion(uTags otherNGS)try :uGenericNGS(otherNGS.getChr(),otherNGS.getStart(),otherNGS.getEnd(), otherNGS.getStrand())
{}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(uTags)");
    e << region_error(*this);
    throw e;
}
uRegion::uRegion(uBasicNGS otherNGS)try :uGenericNGS(otherNGS.getChr(),otherNGS.getStart(),otherNGS.getEnd(), otherNGS.getStrand())
{}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(uBasicNGS)");
    e << region_error(*this);
    throw e;
}

  /** \brief Constructor from parser Token
   * \param uToken Valid Token with data. We assumed the token is valid, so skip some checks here to avoid duplication
   */
uRegion::uRegion(uToken pToken)try :uGenericNGS(pToken){


  	 if (pToken.isParamSet(token_param::DENSITY))
                setDensity(utility::stof(pToken.getParam(token_param::DENSITY)));


}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(uToken)");
    e << region_error(*this);
    throw e;
}



/** \brief Destructor
 */
uRegion::~uRegion()
{}

/** \brief Set the signal value for our region at a specific position.
 *
 * \param i int : Position du signal à changer
 * \param value float : Signal value
 * \return void
 *
 */
void uRegion::setSignal(int i, float value)
{
 try{
    if (signal.size()==0)
        signal.resize(this->getLenght());

/**< Sanity check */
    if (i < this->getLenght())
        signal.at(i)=value;
    else
        throw param_throw()<<string_error("Failling in setSignal, trying to set signal outside of region boundary at position "+utility::to_string(i)+"\n");
    }
    catch(param_throw &e){
        #ifdef DEBUG
               cerr << "Throwing in uTags setSignal" <<endl;
        #endif
        e << region_error(*this);
        throw e;
    }
}

bool uRegion::isEqual(const uRegion & pCompared)const{

    return ((this->getChr()==pCompared.getChr())&&
            (this->getStrand()==pCompared.getStrand())&&
            (this->getStart()==pCompared.getStart())&&
            (this->getEnd()==pCompared.getEnd())&&
            (this->getScoreVector()==pCompared.getScoreVector())&&
            (this->getIdent()==pCompared.getIdent())&&
            (this->getDensity()==pCompared.getDensity())&&
            (this->getCount()==pCompared.getCount())&&
            (this->getSignal()==pCompared.getSignal()));

}

/** \brief Set the signal of an entire region, must be appropriate size
 *
 * \param ourSignal std::vector<float>
 * \return void
 *
 */
void uRegion::setSignal(std::vector<float> ourSignal)
{
    try {
    if ((int)ourSignal.size()<getLenght())
        throw param_throw();
    signal = ourSignal;
    }
    catch(param_throw & e)
    {
        #ifdef DEBUG
               cerr << "Failling in uRegion setSignal(vector), received a vector of size greater then elem lenght of "+utility::to_string(getLenght())+"\n" <<endl;
        #endif
        e<<string_error("Failling in uRegion setSignal(vector), received a vector of size greater then elem lenght of "+utility::to_string(getLenght())+"\n");
        e<<region_error(*this);
        throw e;
    }
}

/** \brief Return ou signal vector for the region, may be empty
 *
 * \return std::vector<float> :: Signal returned
 *
 */
std::vector<float> uRegion::getSignal() const
{
    return signal;
}

/** \brief Output our signal data
 *
 * \param out std::ostream& : Output stream
 * \return void
 *
 */
 //TODO use parser for all write functions
//void uRegion::writeAll(std::ostream& out) const
//{
//    std::vector<float>::iterator iterVec;
//
//    out << getChr()<<"\t" << getStart() << "\t"<< getEnd() << "\t" << getCount();
//
//    for (int i=0; i< getScoreCount() ; i++)
//        out <<  "\t" << getScore(i);
//
//    out << std::endl;
//}

/** \brief Output our signal data with start/end
 *
 * \param out std::ostream& : Output stream
 * \return void
 *
 */
void uRegion::writeSignal(std::ostream& out) const
{
  //  out << getChr() << "\t" << getStart() << "\t" << getEnd() <<  endl;
    for (auto iterVec = signal.begin() ; iterVec!= signal.end(); iterVec++)
    {
        out << *iterVec <<  "\t";
    }
    out << std::endl;
}

/** \brief Output our region, minus signal
 *
 * \param out std::ostream& : Ofstream
 * \return void
 *
 */
//void uRegion::writeRegion(std::ostream& out) const
//{
//    out << getChr()<<"\t" << getStart() << "\t"<< getEnd() << "\t" << getCount() <<   std::endl;
//}

uRegion uRegion::getCopy() const{

        uRegion copyElem =  *this;
        return copyElem;
}

/**< Chrom constructors */
    uRegionChrom::uRegionChrom(uBasicNGSChrom pCopyChrom):uGenericNGSChrom(pCopyChrom.getChr()){
        chromSize=pCopyChrom.getChromSize();
        for (auto itr= pCopyChrom.begin(); itr!=pCopyChrom.end(); itr++  )
            addData(uRegion(*itr));
    }

    uRegionChrom::uRegionChrom(uTagsChrom pCopyChrom):uGenericNGSChrom(pCopyChrom.getChr()){
            chromSize=pCopyChrom.getChromSize();
            for (auto itr= pCopyChrom.begin(); itr!=pCopyChrom.end(); itr++  )
                addData(uRegion(*itr));
    }

    uRegionChrom::uRegionChrom(const uGenericNGSChrom<uRegionChrom,uRegion> & copyCop)
    {
        VecSites=copyCop.returnVecData();
        chr= copyCop.getChr();
        m_isSorted=copyCop.getSortedStatus();
        sortGetStart=copyCop.getStartFunct();
        sortGetEnd=copyCop.getEndFunct();
        m_comptFunc=copyCop.getCompFunct();
        chromSize=copyCop.getChromSize();
    }
    uRegionChrom::uRegionChrom(const uRegionChrom& initFrom){

        VecSites=initFrom.returnVecData();
        chr= initFrom.getChr();
        m_isSorted=initFrom.m_isSorted;
        sortGetStart=initFrom.sortGetStart;
        sortGetEnd=initFrom.sortGetEnd;
        m_comptFunc=initFrom.m_comptFunc;
        chromSize=initFrom.chromSize;

    }
    uRegionChrom& uRegionChrom::operator=(const uRegionChrom& copFrom)
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

    uRegionChrom uRegionChrom::getCopy() const{

        uRegionChrom copyElem =  *this;
        return copyElem;
    }

//TODO Once implement, this needs to use our Transform wrapper
/** \brief Set the density scores for our exp versus another one
    Overlaps exist for uTags, uRegion, uBase
 *
 * \param : The experiment we want to get our density comparison from
 * \param OverlapType const : Type of overlap we use
 * \return void
 *
 */

void uRegionExperiment::measureDensityOverlap(const uRegionExperiment& expToComp, const OverlapType poverlap)
{
    for (auto it =ExpMap.begin(); it!=ExpMap.end(); it++ )
    {
        const uRegionChrom* chromToCompare=expToComp.getpChrom(it->second.getChr());
        (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
    }
}

void uRegionExperiment::measureDensityOverlap(const uTagsExperiment& expToComp, const OverlapType poverlap)
{
    for (auto it =ExpMap.begin(); it!=ExpMap.end(); it++ )
    {
        if (expToComp.getpChrom(it->second.getChr())!=nullptr){
            const uTagsChrom* chromToCompare=expToComp.getpChrom(it->second.getChr());
            (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
            }
    }
}
void uRegionExperiment::measureDensityOverlap(const uBasicNGSExperiment& expToComp, const OverlapType poverlap)
{
    for (auto it =ExpMap.begin(); it!=ExpMap.end(); it++ )
    {
        if (expToComp.getpChrom(it->second.getChr())!=nullptr){
            const uBasicNGSChrom* chromToCompare=expToComp.getpChrom(it->second.getChr());
            (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
            }
    }
}

/** \brief Output all signal data of Exp
 * \param out std::ostream& : Output stream
 * \return void
 *
 */
void uRegionExperiment::writeSignal(std::ostream& out)
{
    applyOnSites(bind2nd(mem_fun_ref(&uRegion::writeSignal), out));
}

/** \brief  Set the density scores for our chrom versus another one
 *             OVerloads exist for uTags, uRegion, uBase
 * \param uGenericNGSChrom& : Chrom to get density comparison from
 * \param OverlapType const : Type of overlap we use
 * \return void
 *
 */

void uRegionChrom::measureDensityOverlap( const uBasicNGSChrom& chromtoComp, const OverlapType pOverlap)
{
    for (auto it =VecSites.begin(); it!=VecSites.end(); it++ )
    {
        it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}
void uRegionChrom::measureDensityOverlap(const uTagsChrom& chromtoComp, const OverlapType pOverlap)
{
    for (auto it =VecSites.begin(); it!=VecSites.end(); it++ )
    {
        it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}
void uRegionChrom::measureDensityOverlap(const uRegionChrom& chromtoComp, const OverlapType pOverlap)
{
    for (auto it =VecSites.begin(); it!=VecSites.end(); it++ )
    {
        it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}

//TODO re-write so it does not required a density vector.
 /** \brief Generate Density signal from another region Chrom
  *         Overloads exist for uRegion, uTags, uBase
  * \param chromToComp const uRegionChrom&
  * \return void
  *
  */
 void uRegionChrom::generateSignal(const uRegionChrom & chromToComp)
{
    vector<long int> densityValues;
    string trace;
try
    {
        densityValues.resize(chromToComp.getChromSize());
        trace += ("Starting apply on All Sites (uRegion) \n");
        vector<uRegion> failRegionVec;
        chromToComp.applyOnAllSites([&] (const uRegion & Elem)
        {
            try
            {
                /**< Ignore if we map over the reference */
                if (Elem.getEnd() <(int)densityValues.size())
                    for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
                    {
                        densityValues.at(i)+=Elem.getCount();
                    }
                else
                {
                    /**< Keep log of invalid Regions */
                    failRegionVec.push_back(Elem);
                }
            }
            catch (...)
            {
                throw elem_throw() << region_error(Elem);
            }
        }
                               );
        trace += ("Starting apply on All Sites (Uregion) \n");
        this->applyOnAllSites( [&] (uRegion& Elem)
        {
            try
            {
                vector<float> signalVector;
                signalVector.resize(Elem.getLenght());
                /**< Ignore if over the reference */
                if (Elem.getEnd() <(int)densityValues.size())
                {
                    for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
                    {
                        signalVector.at(Elem.getEnd()-i)=densityValues.at(i);
                    }
                    Elem.setSignal(signalVector);
                }
                else
                {
                    failRegionVec.push_back(Elem);
                    cerr << "Not generating signal for region, as it exceeds reference size of"<< (int)densityValues.size() <<endl;

                }
            }
            catch(...)
            {
                throw elem_throw() << region_error(Elem);
            }
        } );
        if (failRegionVec.size())
        {
            skipped_elem_throw e;
            if (failRegionVec.size())
                e <<skipped_regions(failRegionVec);
            e << string_error(trace);
            throw e;
        }
    }
    catch(skipped_elem_throw & e)
    {
        throw e;
    }
    catch(elem_throw & e)
    {
        e << string_error("Throwing in uRegionChrom::generateSignal(uTagsExperiment& expToComp), from elem_throw error \n Reference is size at "+(utility::to_string(densityValues.size()))+
                          "\n"+
                          trace);
        throw e;
    }
    catch(std::exception & e)
    {
        throw ugene_exception_base()
                << string_error("Throwing in uRegionChrom::generateSignal(uTagsExperiment& expToComp), from out_of_range error \n"+trace);
    }

}

void uRegionChrom::generateSignal(const uTagsChrom & chromToComp)
{
    vector<long int> densityValues;
    string trace;
try
    {
        densityValues.resize(chromToComp.getChromSize());
        trace += ("Starting apply on All Sites (uRegion) \n");
        vector<uTags> failTagVec;
        chromToComp.applyOnAllSites([&] (const uTags & Elem)
        {
            try
            {
                /**< Ignore if we map over the reference */
                if (Elem.getEnd() <(int)densityValues.size())
                    for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
                    {
                        densityValues.at(i)+=1;
                    }
                else
                {
                    /**< Keep log of invalid Regions */
                    failTagVec.push_back(Elem);
                }
            }
            catch (...)
            {
                throw elem_throw() << tag_error(Elem);
            }
        }
                               );
         vector<uRegion> failRegionVec;
        trace += ("Starting apply on All Sites (Uregion) \n");
        this->applyOnAllSites( [&] (uRegion& Elem)
        {
            try
            {
                vector<float> signalVector;
                signalVector.resize(Elem.getLenght());
                /**< Ignore if over the reference */
                if (Elem.getEnd() <(int)densityValues.size())
                {
                    for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
                    {
                        signalVector.at(Elem.getEnd()-i)=densityValues.at(i);
                    }
                    Elem.setSignal(signalVector);
                }
                else
                {
                    failRegionVec.push_back(Elem);
                    cerr << "Not generating signal for region, as it exceeds reference size of"<< (int)densityValues.size() <<endl;

                }
            }
            catch(...)
            {
                throw elem_throw() << region_error(Elem);
            }
        } );
        if (( failRegionVec.size())||(failTagVec.size()))
        {
            skipped_elem_throw e;
            if (failRegionVec.size())
                e <<skipped_regions(failRegionVec);
            if (failTagVec.size())
                e <<skipped_tags(failTagVec);
            e << string_error(trace);
            throw e;
        }
    }
    catch(skipped_elem_throw & e)
    {
        throw e;
    }
    catch(elem_throw & e)
    {
        e << string_error("Throwing in uRegionChrom::generateSignal(uTagsExperiment& expToComp), from elem_throw error \n Reference is size at "+(utility::to_string(densityValues.size()))+
                          "\n"+
                          trace);
        throw e;
    }
    catch(std::exception & e)
    {
        throw ugene_exception_base()
                << string_error("Throwing in uRegionChrom::generateSignal(uTagsExperiment& expToComp), from out_of_range error \n"+trace);
    }

}

void uRegionChrom::generateSignal(const uBasicNGSChrom& chrToComp)
{

//TODO redo error management
    vector<long int> densityValues;
    string trace;

    try
    {
        densityValues.resize(chrToComp.getChromSize());
        trace += ("Starting apply on All Sites (Utags) \n");
        vector<uBasicNGS> failBasicvec;
        (chrToComp).applyOnAllSites([&] (const uBasicNGS & Elem)
        {
            try
            {
                /**< Ignore if we map over the reference */
                if (Elem.getEnd() <(int)densityValues.size())
                    for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
                    {
                        densityValues.at(i)++;
                    }
                else
                {
                    /**< Keep log of invalid tags */
                    failBasicvec.push_back(Elem);
                }
            }
            catch (...)
            {
                throw elem_throw() << basic_error(Elem);
            }
        }
    );
        trace += ("Starting apply on All Sites (Uregion) \n");
        vector<uRegion> failRegionVec;
        this->applyOnAllSites( [&] (uRegion& Elem)
        {
            try
            {
                vector<float> signalVector;
                signalVector.resize(Elem.getLenght());
                /**< Ignore if over the reference */
                if (Elem.getEnd() <(int)densityValues.size())
                {
                    for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
                    {
                        signalVector.at(Elem.getEnd()-i)=densityValues.at(i);
                    }
                    Elem.setSignal(signalVector);
                }
                else
                {
                    failRegionVec.push_back(Elem);
                }
            }
            catch(...)
            {
                throw elem_throw() << region_error(Elem);
            }
        } );

        if (failBasicvec.size() || (failRegionVec.size()))
        {
            skipped_elem_throw e;
            if (failBasicvec.size())
                e <<skipped_Basic(failBasicvec);
            if (failRegionVec.size())
                e <<skipped_regions(failRegionVec);
            e << string_error(trace);
            throw e;
        }
    }
    catch(skipped_elem_throw & e)
    {
        throw e;
    }
    catch(elem_throw & e)
    {
        e << string_error("Throwing in uRegionChrom::generateSignal(uTagsExperiment& expToComp), from elem_throw error \n Reference is size at "+
                         (utility::to_string((int)densityValues.size()))+
                          "\n"+
                          trace);
        throw e;
    }
    catch(std::exception & e)
    {
        throw ugene_exception_base()
                << string_error("Throwing in uRegionChrom::generateSignal(uTagsExperiment& expToComp), from out_of_range error \n"+trace);
    }
}


/** \brief Get number of tags overlapping each position of each elements of each chrom
    Overloads exist for uTags, uRegion, uBase
 *
 * \param out std::ostream&
 * \return void
 *
 */
void uRegionExperiment::generateSignal(const uTagsExperiment& expToComp)
{
    vector<uTags> skipTagvec;
    vector<uRegion> skipRegioNvec;

    for(auto& chrom : ExpMap)
    {
       try {
            const uTagsChrom* chromToCompare=expToComp.getpChrom(chrom.second.getChr());
            chrom.second.generateSignal(*chromToCompare);
        }
        /**<  Managed invalid tags. Eventually, this should offer HARD/SOFT/SILENT Options */
        catch(skipped_elem_throw & e)
        {
            if (vector<uTags> const * vecU =boost::get_error_info<skipped_tags>(e) )
                {
                    skipTagvec.insert( skipTagvec.end(), vecU->begin(), vecU->end() );
                }
               if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
               {
                    skipRegioNvec.insert( skipRegioNvec.end(), vecR->begin(), vecR->end() );
               }
        }
    }
/**< If we caught any errors, we log them */
   if ( (skipTagvec.size()>0) || (skipRegioNvec.size()>0))
        {
            skipped_elem_throw e;
            if (skipTagvec.size())
                e <<skipped_tags(skipTagvec);
            if (skipRegioNvec.size())
                e <<skipped_regions(skipRegioNvec);
            throw e;
        }
}
 void uRegionExperiment::generateSignal(const uRegionExperiment & expToComp){

 vector<uTags> skipTagvec;
 vector<uRegion> skipRegioNvec;

for(auto& chrom : ExpMap)
    {
       try {
            const uRegionChrom* chromToCompare=expToComp.getpChrom(chrom.second.getChr());
            chrom.second.generateSignal(*chromToCompare);
        }
        catch(skipped_elem_throw & e)
        {
            if (vector<uTags> const * vecU =boost::get_error_info<skipped_tags>(e) )
                {
                  skipTagvec.insert( skipTagvec.end(), vecU->begin(), vecU->end() );
                }
               if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
               {
                    skipRegioNvec.insert( skipRegioNvec.end(), vecR->begin(), vecR->end() );
               }
        }
    }
/**< Log errors */
 if ( (skipTagvec.size()>0) || (skipRegioNvec.size()>0))
        {
            skipped_elem_throw e;
            if (skipTagvec.size())
                e <<skipped_tags(skipTagvec);
            if (skipRegioNvec.size())
                e <<skipped_regions(skipRegioNvec);
            throw e;
        }
 }
 void uRegionExperiment::generateSignal(const uBasicNGSExperiment & expToComp){
 vector<uBasicNGS> skipBasic;
 vector<uRegion> skipRegioNvec;
for(auto& chrom : ExpMap)
    {
       try {
            const uBasicNGSChrom* chromToCompare=expToComp.getpChrom(chrom.second.getChr());
            chrom.second.generateSignal(*chromToCompare);
        }
        catch(skipped_elem_throw & e)
        {
            if (vector<uBasicNGS> const * vecU =boost::get_error_info<skipped_Basic>(e) )
                {
                  skipBasic.insert( skipBasic.end(), vecU->begin(), vecU->end() );
                }
               if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
               {
                    skipRegioNvec.insert( skipRegioNvec.end(), vecR->begin(), vecR->end() );
               }
        }
    }
/**< Log skipBasic */
 if ( (skipBasic.size()>0) || (skipRegioNvec.size()>0))
        {
            skipped_elem_throw e;
            if (skipBasic.size())
                e <<skipped_Basic(skipBasic);
            if (skipRegioNvec.size())
                e <<skipped_regions(skipRegioNvec);
            throw e;
        }
 }


/** \brief Write all chromosome data to a stream
 * \param std::ostream&: Output stream to use.
 */
//void uRegionChrom::writeAll(std::ostream& out)
//{
//    for(auto& it : VecSites)
//    {
//        it.writeAll(out);
//    }
//}

/** \brief Write all the experiment to a stream
 * \param std::ostream&: Output stream to use.
 */
//void uRegionExperiment::writeAll(std::ostream& out)
//{
//    for(auto& chrom : ExpMap)
//    {
//        chrom.second.writeAll(out);
//    }
//}

/** \brief Write the experiment's density data as tab separated values
 * \param std::ostream&: Output stream to use.
 */
//void uRegionExperiment::writeDensityAsTab(std::ostream& out)
//{
//    for(auto& chrom : ExpMap)
//    {
//        chrom.second.writeDensityAsTab(out);
//    }
//}

/** \brief Write the experiment's density data as tab separated values
 * \param std::ostream&: Output stream to use.
 */
//void uRegionChrom::writeDensityAsTab(std::ostream& out)
//{
//    for(auto& it : VecSites)
//    {
//        it.writeRegion(out);
//    }
//}


uRegionExperiment uRegionExperiment::getCopy() const{

    uRegionExperiment copyElem =  *this;
    return copyElem;
}


//TODO Test our loader for Wig
//TODO Move to parser?
/** \brief Load our density data from a wig file
 *
 * \param std::ifstream valid Input stream
 * \return void
 *
 */
void uRegionExperiment::loadFromWig(std::ifstream & inputStream)
{
try {
string curChr, curLine;
int curSpan,curStart, curStep;
enum class wigType
{
    VARIABLE_STEP, FIXED_STEP, NONE
};

/**<  Validate header */
getline(inputStream, curLine);
if (curLine.find("track")==string::npos)
    throw format_parsing_error() <<string_error("No track definition line in wiggle file\n");

/**< At start, we do not know what to process */
wigType curWig = wigType::NONE;

/**<Load line  */
    while(!getline(inputStream, curLine).eof()){
{
        utility::Tokenizer data(curLine);
        data.NextToken();
        string firstToken= data.GetToken();

        /**< Is it a track definition line? Format forces the following two tags */
        if((firstToken=="variableStep")||(firstToken=="variableStep"))
        {
            /**< Parse */
            if (firstToken=="variableStep")
            {
                curWig=wigType::VARIABLE_STEP;
                /**< Missing chrom tag */
                if (!(data.NextToken()))
                        throw 17;

                curChr=data.GetToken();
                /**< If invalid header */
                if (curChr.find("chrom=")==string::npos)
                    throw 17;

                curChr=curChr.substr(curChr.find("chrom="));

                /**< Optional Span parameter */
                 curSpan=0;
                 if(data.NextToken()){
                    string span;
                    span = data.GetToken();
                    /**< If good, yay, if not, fail again*/
                    if (span.find("span=")==string::npos)
                        throw 18;

                     span=span.substr(span.find("span="));
                     curSpan=utility::stoi(span);
                 }
            }/**< Fixed Step */
            else
            {
                 curWig=wigType::FIXED_STEP;
                /**< Chrom */
                 if (!(data.NextToken()))
                        throw 19;
                  curChr=data.GetToken();
                  if (curChr.find("chrom=")==string::npos)
                        throw 19;
                  curChr=curChr.substr(curChr.find("chrom="));
                    /**< Start */
                  if (!(data.NextToken()))
                        throw 20;

                  string start= data.GetToken();
                  if (start.find("start=")==string::npos)
                        throw 20;

                   start=start.substr(start.find("start="));
                    curStart=utility::stoi(start);
                    /**< Step */
                    if (!(data.NextToken()))
                        throw 21;

                  string step= data.GetToken();
                  if (step.find("step=")==string::npos)
                        throw 21;

                   step=step.substr(start.find("step="));
                   curStep=utility::stoi(step);

                    /**<Span  */

                    /**< Optional Span parameter */
                    curSpan=0;
                 if(data.NextToken()){
                    string span;
                    span = data.GetToken();
                    /**< If good, yay, if not, fail again*/
                    if (span.find("span=")==string::npos)
                        throw 22;

                     span=span.substr(span.find("span="));
                     curSpan=utility::stoi(span);
                 }
            }
        }
        else
        {
        /**< Process the line */
         uRegion curReg;
         switch (curWig)
            {
                /**< Crash */
            case wigType::NONE:
                cerr << "First line in a wig file must be track declaration, failling" <<endl;
                throw 23;
                break;
            case wigType::FIXED_STEP:
                curReg.setChr(curChr);
                curReg.setStartEnd(curStart,curStart+curSpan);
                curReg.setCount(utility::stoi(firstToken));
                curStart+=curStep;
                break;
            case  wigType::VARIABLE_STEP:

                 data.NextToken();
                 int tagcount=utility::stoi(data.GetToken());
                 curReg.setChr(curChr);
                 int start=utility::stoi(firstToken);
                 curReg.setStartEnd(start,start+curSpan);
                 curReg.setCount(tagcount);
                break;

            }
            addData(curReg);
        }
}

    inferChrSize();
    }

}
    catch(format_parsing_error & e )
    {
        cerr << "Catching while trying to load wiggle file, format error" <<endl;
        throw;
    }
    catch(...)
    {
        cerr << "Catching while trying to load wiggle file, passing" <<endl;
        throw;
    }

}
} // End of namespace NGS
