#include "uFormats.h"
#include "uRegion.h"
#include "uTags.h"

using namespace std;

/** \brief Default constructor
 */
uRegion::uRegion()
{

}

/** \brief Region constructor, with start, end and chr
 *
 * \param ourchr std::string  : Chrom to set
 * \param ourstart int        : Start position
 * \param ourend int          : End position
 *
 */
uRegion::uRegion(std::string ourchr, int ourstart, int ourend)
{
    setChr(ourchr);
    setStartEnd(ourstart,ourend);
}


/** \brief Constructor from parent singular
 *
 * \param uGenericNGS
 *
 */
uRegion::uRegion(uGenericNGS otherNGS)
{
    setChr(otherNGS.getChr());
    setStartEnd(otherNGS.getStart(),otherNGS.getEnd());
}


uRegion::~uRegion()
{
    //dtor
}

void uRegion::setScore(float p_score, int p_Pos)
{
    if (p_Pos>= (score.size()))
        score.resize(p_Pos+1);
    score.at(p_Pos)=p_score;
}

float uRegion::getScore(int p_Pos) const
{
    if (p_Pos>= (score.size()))
       return 0;
    return score.at(p_Pos);
}


/** \brief Set the signal value for our region at a specific position.
 *
 * \param i int : Position du signal à changer
 * \param value float : Signal value
 * \return void
 *
 */
void uRegion::setSignal(int i, float value)
{
    if (signal.size()==0)
        signal.resize(this->getLenght());

//Sanity check
//TODO Explicit error message
    if (i < this->getLenght())
        signal.at(i)=value;
}

/** \brief Set the signal of an entire region
 *
 * \param ourSignal std::vector<float>
 * \return void
 *
 */
void uRegion::setSignal(std::vector<float> ourSignal)
{
    signal= ourSignal;
}

/** \brief Return ou signal vector for the region
 *
 * \return std::vector<float> :: Signal returned
 *
 */
std::vector<float> uRegion::getSignal()
{

    return signal;
}

/** \brief Output our signal data
 *
 * \param out std::ostream& : Output stream
 * \return void
 *
 */
void uRegion::writeAll(std::ostream& out) const
{

    std::vector<float>::iterator iterVec;

    out << getChr()<<"\t" << getStart() << "\t"<< getEnd() << "\t" << getCount();

    for (int i=0; i< getScoreCount() ; i++)
        out <<  "\t" << getScore(i);

    out << std::endl;
    //for (iterVec = signal.begin() ;iterVec!= signal.end(); iterVec++){
    //    out << *iterVec <<  " ";
    //    }
}


/** \brief Output our signal data with start/end
 *
 * \param out std::ostream& : Output stream
 * \return void
 *
 */
void uRegion::writeSignal(std::ostream& out) const
{
    out << getChr() << "\t" << getStart() << "\t" << getEnd() <<  endl;
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
void uRegion::writeRegion(std::ostream& out) const
{
    out << getChr()<<"\t" << getStart() << "\t"<< getEnd() << "\t" << getCount() <<   std::endl;
}


//ToDO Once implement, this needs to use our Transform wrapper
/** \brief Set the density scores for our exp versus another one
 *
 * \param uGenericNGSExperiment : The experiment we want to get our density comparison from
 * \param OverlapType const : Type of overlap we use
 * \return void
 *
 */
void uRegionExperiment::measureDensityOverlap(uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>,uGenericNGS>& expToComp, const OverlapType poverlap)
{
    for (auto it =ExpMap.begin(); it!=ExpMap.end(); it++ )
    {
        uGenericNGSChrom<uGenericNGS>* chromToCompare=expToComp.getpChrom(it->second.getChr());
        (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
    }
}

void uRegionExperiment::measureDensityOverlap(uTagsExperiment& expToComp, const OverlapType poverlap)
{
    for (auto it =ExpMap.begin(); it!=ExpMap.end(); it++ )
    {
        uTagsChrom* chromToCompare=expToComp.getpChrom(it->second.getChr());
        (it)->second.measureDensityOverlap(*chromToCompare,poverlap);
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
 *
 * \param uGenericNGSChrom& : Chrom to get density comparison from
 * \param OverlapType const : Type of overlap we use
 * \return void
 *
 */
void uRegionChrom::measureDensityOverlap(uGenericNGSChrom<uGenericNGS>& chromtoComp, const OverlapType pOverlap)
{
    for (auto it =VecSites.begin(); it!=VecSites.end(); it++ )
    {
        it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}

void uRegionChrom::measureDensityOverlap(uTagsChrom& chromtoComp, const OverlapType pOverlap)
{
    for (auto it =VecSites.begin(); it!=VecSites.end(); it++ )
    {
        it->setCount((chromtoComp.getSubsetCount(it->getStart(), it->getEnd(),pOverlap)));
    }
}



/** \brief  Set the density scores for our chrom versus another one
 *
 * \param uGenericNGSExperiment& : Exp to get density comparison from
 * \param OverlapType const : Type of overlap we use
 * \return void
 *
 */
void uRegionChrom::measureDensityOverlap(uGenericNGSExperiment< uGenericNGSChrom<uGenericNGS>,uGenericNGS>& expToCompare, const OverlapType pOverlap)
{
    uGenericNGSChrom<uGenericNGS>* chromToCompare;
    chromToCompare=  expToCompare.getpChrom(this->getChr());
    measureDensityOverlap(*chromToCompare, pOverlap);
}



/** \brief Get number of tags overlapping each position of each elements. note that we ignore tags over the reference.
 *
 * \param expToComp uTagsExperiment&
 * \return void
 *
 */
void uRegionChrom::generateSignal(const uRegionExperiment & expToComp){

    const uRegionChrom *  pChrom;
    try {
        pChrom=expToComp.getpChrom(getChr());
        this->generateSignal(*pChrom);
    }
    catch(...){
        throw;
    }
}

//TODO re-write so it does not required a density vector.
 /** \brief Generate Density signal from another region Chrom
  *
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
        vector<float> densityVector;

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
                throw elem_throw() << tag_error(Elem);
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
                 //  Elem.debugElem();
                 //  utility::pause_input();
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
        e << string_error("Throwing in uRegionChrom::generateSignal(uTagsExperiment& expToComp), from elem_throw error \n Reference is size at "+
                         (utility::numberToString((int)densityValues.size()))+
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


/** \brief Get number of tags overlapping each position of each elements. note that we ignore tags over the reference.
 *
 * \param expToComp uTagsExperiment&
 * \return void
 *
 */
void uRegionChrom::generateSignal(uTagsExperiment& expToComp)
{

//TODO redo error management
     vector<long int> densityValues;
    string trace;
    try
    {

        auto pTag =expToComp.getpChrom(this->getChr());
        densityValues.resize(pTag->getChromSize());

        trace += ("Starting apply on All Sites (Utags) \n");
        vector<uTags> failTagvec;
        (pTag)->applyOnAllSites([&] (uTags Elem)
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
                   // cerr << "Ignore following tags as maps over reference of size "<< (int)densityValues.size() <<endl;
                    /**< Keep log of invalid tags */
                    failTagvec.push_back(Elem);
                   // Elem.debugElem();
                }
            }
            catch (...)
            {
                throw elem_throw() << tag_error(Elem);
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
                    cerr << "Not generating signal for region, as it exceeds reference size of"<< (int)densityValues.size() <<endl;
                 //  Elem.debugElem();
                 //  utility::pause_input();
                }
            }
            catch(...)
            {
                throw elem_throw() << region_error(Elem);
            }
        } );

        if (failTagvec.size() || (failRegionVec.size()))
        {
            skipped_elem_throw e;
            if (failTagvec.size())
                e <<skipped_tags(failTagvec);
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
                         (utility::numberToString((int)densityValues.size()))+
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
 *
 * \param out std::ostream&
 * \return void
 *
 */
void uRegionExperiment::generateSignal(uTagsExperiment& expToComp)
{
    vector<uTags> skipTagvec;
    vector<uRegion> skipRegioNvec;

    for(auto& chrom : ExpMap)
    {
       try {
        chrom.second.generateSignal(expToComp);
        }
        /**<  Managed invalid tags. Eventually, this should offer HARD/SOFT/SILENT Options */
        catch(skipped_elem_throw & e)
        {
            // cerr << "Catching at Reg level, was testing " << chrom.second.getChr() <<endl;
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
/**< Generate Density overlap from another density Map */
 void uRegionExperiment::generateSignal(const uRegionExperiment & expToComp){

 vector<uRegion> skipTagvec;
 vector<uRegion> skipRegioNvec;

for(auto& chrom : ExpMap)
    {
       try {
        chrom.second.generateSignal(expToComp);
        }
        catch(skipped_elem_throw & e)
        {
            // cerr << "Catching at Reg level, was testing " << chrom.second.getChr() <<endl;
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
                e <<skipped_regions(skipTagvec);
            if (skipRegioNvec.size())
                e <<skipped_regions(skipRegioNvec);
            throw e;
        }

 }


void uRegionChrom::writeAll(std::ostream& out)
{
    for(auto& it : VecSites)
    {
        it.writeAll(out);
    }

}

/**<  */

void uRegionExperiment::writeAll(std::ostream& out)
{
    for(auto& chrom : ExpMap)
    {
        chrom.second.writeAll(out);
    }
}

/** @brief writeDensityAsTab
  *
  * @todo: document this function
  */
void uRegionExperiment::writeDensityAsTab(std::ostream& out)
{
    for(auto& chrom : ExpMap)
    {
        chrom.second.writeDensityAsTab(out);
    }
}

/** @brief writeDensityAsTab
  *
  * @todo: document this function
  */
void uRegionChrom::writeDensityAsTab(std::ostream& out)
{
    for(auto& it : VecSites)
    {
        it.writeRegion(out);
    }


}
//TODO Test our loader for Wig
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
                     curSpan=utility::stringToInt(span);
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
                    curStart=utility::stringToInt(start);
                    /**< Step */
                    if (!(data.NextToken()))
                        throw 21;

                  string step= data.GetToken();
                  if (step.find("step=")==string::npos)
                        throw 21;

                   step=step.substr(start.find("step="));
                   curStep=utility::stringToInt(step);

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
                     curSpan=utility::stringToInt(span);
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
                curReg.setCount(utility::stringToInt(firstToken));
                curStart+=curStep;
                break;
            case  wigType::VARIABLE_STEP:
                 data.NextToken();
                 int tagcount=utility::stringToInt(data.GetToken());
                 curReg.setChr(curChr);
                 int start=utility::stringToInt(firstToken);
                 curReg.setStartEnd(start,start+curSpan);
                 curReg.setCount(tagcount);
                break;

            }
            addSite(curReg);
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


