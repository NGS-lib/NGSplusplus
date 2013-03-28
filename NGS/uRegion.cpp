
// ***************************************************************************
// uRegion.cpp (c) 2013
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
// **************************************************************************


#include "uFormats.h"
#include "uRegion.h"
#include "uTags.h"
#include "uBasicNGS.h"
#include "IO/uToken.h"
#include "uGene.h"

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
{
}
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
{
 setScoreVector(otherNGS.getScoreVector());
}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(uTags)");
    e << region_error(*this);
    throw e;
}
uRegion::uRegion(uBasicNGS otherNGS)try :uGenericNGS(otherNGS.getChr(),otherNGS.getStart(),otherNGS.getEnd(), otherNGS.getStrand())
{
 setScoreVector(otherNGS.getScoreVector());
}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(uBasicNGS)");
    e << region_error(*this);
    throw e;
}

uRegion::uRegion(uGene otherNGS)try :uGenericNGS(otherNGS.getChr(),otherNGS.getStart(),otherNGS.getEnd(), otherNGS.getStrand())
{
 setScoreVector(otherNGS.getScoreVector());
}
catch(construct_elem_throw &e)
{
    addStringError(e,"Throwing in uRegion(uGene)");
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
    if ((int)ourSignal.size()!=getLenght())
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


/** \brief Output our signal data, with seperator as requied
 *
 * \param out std::ostream& : Output stream
 * \return void
 *
 */
void uRegion::writeSignal(std::ostream& out, char pSep) const
{
  //  out << getChr() << "\t" << getStart() << "\t" << getEnd() <<  endl;
  if (signal.size()){
    for (auto iterVec = signal.begin() ; iterVec!= signal.end(); iterVec++)
    {
        out << *iterVec << pSep;
    }
    out << std::endl;
}
}


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


    uRegionChrom::uRegionChrom(uGeneChrom pCopyChrom):uGenericNGSChrom(pCopyChrom.getChr()){
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

/** \brief Output all signal data of Exp using pSep as seperator
 * \param out std::ostream& : Output stream
 * \return void
 *
 */
void uRegionChrom::writeSignal(std::ostream& out, char pSep)
{
  //  uRegion hihi;
    auto applyFunct=(std::bind(&uRegion::writeSignal,std::placeholders::_1,ref(out),pSep) );
    applyOnAllSites(applyFunct);
    //applyOnAllSites(bind2nd(mem_fun_ref(&uRegion::writeSignal), out));
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
     //   it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
      it->setCount((chromtoComp.getOverlappingCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}
void uRegionChrom::measureDensityOverlap(const uTagsChrom& chromtoComp, const OverlapType pOverlap)
{
    for (auto it =VecSites.begin(); it!=VecSites.end(); it++ )
    {
     //   it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
      it->setCount((chromtoComp.getOverlappingCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}

void uRegionChrom::measureDensityOverlap(const uGeneChrom& chromtoComp, const OverlapType pOverlap)
{
    for (auto it =VecSites.begin(); it!=VecSites.end(); it++ )
    {
    //    it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
     it->setCount((chromtoComp.getOverlappingCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}

void uRegionChrom::measureDensityOverlap(const uRegionChrom& chromtoComp, const OverlapType pOverlap)
{
    for (auto it =VecSites.begin(); it!=VecSites.end(); it++ )
    {
       // it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
         it->setCount((chromtoComp.getOverlappingCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}

 /** \brief Generate Density signal from another region Chrom. This signal is stored internaly in the signal value of each element.
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
                if (Elem.getEnd() <=(int)densityValues.size())
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
                if (Elem.getEnd() <=(int)densityValues.size())
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
                    #ifdef DEBUG
      cerr << "Not generating signal for region, as it exceeds reference size of"<< (int)densityValues.size() <<endl;
#endif


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
                if (Elem.getEnd() <=(int)densityValues.size())
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
                if (Elem.getEnd() <=(int)densityValues.size())
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
                    #ifdef DEBUG
      cerr << "Not generating signal for region, as it exceeds reference size of"<< (int)densityValues.size() <<endl;
#endif
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
                if (Elem.getEnd() <=(int)densityValues.size())
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
                if (Elem.getEnd() <=(int)densityValues.size())
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


void uRegionChrom::generateSignal(const uGeneChrom& chrToComp)
{
//TODO redo error management
    vector<long int> densityValues;
    string trace;
    try
    {
        densityValues.resize(chrToComp.getChromSize());
        trace += ("Starting apply on All Sites (Utags) \n");
        vector<uGene> failChromvec;
        (chrToComp).applyOnAllSites([&] (const uBasicNGS & Elem)
        {
            try
            {
                /**< Ignore if we map over the reference */
                if (Elem.getEnd() <=(int)densityValues.size())
                    for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
                    {
                        densityValues.at(i)++;
                    }
                else
                {
                    /**< Keep log of invalid tags */
                    failChromvec.push_back(Elem);
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
                if (Elem.getEnd() <=(int)densityValues.size())
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

        if (failChromvec.size() || (failRegionVec.size()))
        {
            skipped_elem_throw e;
            if (failChromvec.size())
                e <<skipped_genes(failChromvec);
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
        if (expToComp.isChrom(it->second.getChr())){
            const uRegionChrom* chromToCompare=expToComp.getpChrom(it->second.getChr());
            (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
        }
    }
}

void uRegionExperiment::measureDensityOverlap(const uTagsExperiment& expToComp, const OverlapType poverlap)
{
    for (auto it =ExpMap.begin(); it!=ExpMap.end(); it++ )
    {
        if (expToComp.isChrom(it->second.getChr())){
            const uTagsChrom* chromToCompare=expToComp.getpChrom(it->second.getChr());
            (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
            }
    }
}
void uRegionExperiment::measureDensityOverlap(const uBasicNGSExperiment& expToComp, const OverlapType poverlap)
{
    for (auto it =ExpMap.begin(); it!=ExpMap.end(); it++ )
    {
        if (expToComp.isChrom(it->second.getChr())){
            const uBasicNGSChrom* chromToCompare=expToComp.getpChrom(it->second.getChr());
            (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
            }
    }
}

void uRegionExperiment::measureDensityOverlap(const uGeneExperiment& expToComp, const OverlapType poverlap)
{
    for (auto it =ExpMap.begin(); it!=ExpMap.end(); it++ )
    {
        if (expToComp.isChrom(it->second.getChr())){
            const uGeneChrom* chromToCompare=expToComp.getpChrom(it->second.getChr());
            (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
            }
    }
}

/** \brief Output all signal data of Exp
 * \param out std::ostream& : Output stream
 * \return void
 *
 */
void uRegionExperiment::writeSignal(std::ostream& out, char pSep)
{
    auto applyFunct=(std::bind(&uRegion::writeSignal,std::placeholders::_1,ref(out),pSep) );
    applyOnSites(applyFunct);
    //applyOnSites(bind2nd(mem_fun_ref(&uRegion::writeSignal), out));
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
           if (expToComp.isChrom(chrom.second.getChr())){
            const uTagsChrom* chromToCompare=expToComp.getpChrom(chrom.second.getChr());
            chrom.second.generateSignal(*chromToCompare);
            }
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
// void uRegionExperiment::generateSignal(const uTagsExperiment & expToComp){
//
// vector<uTags> skipTagvec;
// vector<uRegion> skipRegioNvec;
//
//for(auto& chrom : ExpMap)
//    {
//       try {
//            const uRegionChrom* chromToCompare=expToComp.getpChrom(chrom.second.getChr());
//            chrom.second.generateSignal(*chromToCompare);
//        }
//        catch(skipped_elem_throw & e)
//        {
//            if (vector<uTags> const * vecU =boost::get_error_info<skipped_tags>(e) )
//                {
//                  skipTagvec.insert( skipTagvec.end(), vecU->begin(), vecU->end() );
//                }
//               if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
//               {
//                    skipRegioNvec.insert( skipRegioNvec.end(), vecR->begin(), vecR->end() );
//               }
//        }
//    }
///**< Log errors */
// if ( (skipTagvec.size()>0) || (skipRegioNvec.size()>0))
//        {
//            skipped_elem_throw e;
//            if (skipTagvec.size())
//                e <<skipped_tags(skipTagvec);
//            if (skipRegioNvec.size())
//                e <<skipped_regions(skipRegioNvec);
//            throw e;
//        }
// }


 void uRegionExperiment::generateSignal(const uGeneExperiment& expToComp)
{
    vector<uGene> skipGeneVec;
    vector<uRegion> skipRegioNvec;

    for(auto& chrom : ExpMap)
    {
       try {
           if (expToComp.isChrom(chrom.second.getChr())){
            const uGeneChrom* chromToCompare=expToComp.getpChrom(chrom.second.getChr());
            chrom.second.generateSignal(*chromToCompare);
           }
        }
        /**<  Managed invalid tags. Eventually, this should offer HARD/SOFT/SILENT Options */
        catch(skipped_elem_throw & e)
        {
            if (vector<uGene> const * vecU =boost::get_error_info<skipped_genes>(e) )
                {
                    skipGeneVec.insert( skipGeneVec.end(), vecU->begin(), vecU->end() );
                }
               if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
               {
                    skipRegioNvec.insert( skipRegioNvec.end(), vecR->begin(), vecR->end() );
               }
        }
    }
/**< If we caught any errors, we log them */
   if ( (skipGeneVec.size()>0) || (skipRegioNvec.size()>0))
        {
            skipped_elem_throw e;
            if (skipGeneVec.size())
                e <<skipped_genes(skipGeneVec);
            if (skipRegioNvec.size())
                e <<skipped_regions(skipRegioNvec);
            throw e;
        }
}
 void uRegionExperiment::generateSignal(const uRegionExperiment & expToComp){
try {
 vector<uRegion> skipRegionVecOne;
 vector<uRegion> skipRegioNvec;

for(auto& chrom : ExpMap)
    {
       try {
           if (expToComp.isChrom(chrom.second.getChr())){
            const uRegionChrom* chromToCompare=expToComp.getpChrom(chrom.second.getChr());
            chrom.second.generateSignal(*chromToCompare);
            }
        }
        catch(skipped_elem_throw & e)
        {
            if (vector<uRegion> const * vecU =boost::get_error_info<skipped_regions>(e) )
                {
                  skipRegionVecOne.insert( skipRegionVecOne.end(), vecU->begin(), vecU->end() );
                }
               if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
               {
                    skipRegioNvec.insert( skipRegioNvec.end(), vecR->begin(), vecR->end() );
               }
        }
    }
/**< Log errors */
 if ( (skipRegionVecOne.size()>0) || (skipRegioNvec.size()>0))
        {
            skipped_elem_throw e;
            skipRegionVecOne.insert( skipRegionVecOne.end(), skipRegioNvec.begin(), skipRegioNvec.end() );
            addStringError(e, "Some regions where skipped in uRegionExperiment generateSignal");
            if (skipRegionVecOne.size())
                e <<skipped_regions(skipRegionVecOne);

            throw e;
        }
 }
catch(...)
    {
    throw;
    }

}

 void uRegionExperiment::generateSignal(const uBasicNGSExperiment & expToComp){
 vector<uBasicNGS> skipBasic;
 vector<uRegion> skipRegioNvec;
for(auto& chrom : ExpMap)
    {
       try {
           if (expToComp.isChrom(chrom.second.getChr())){
            const uBasicNGSChrom* chromToCompare=expToComp.getpChrom(chrom.second.getChr());
            chrom.second.generateSignal(*chromToCompare);

            }
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
uRegionExperiment uRegionExperiment::getCopy() const{

    uRegionExperiment copyElem =  *this;
    return copyElem;
}

} // End of namespace NGS
