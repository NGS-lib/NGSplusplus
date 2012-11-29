#include <stdio.h>
#include <string.h>
#include <sstream>
#include "uTags.h"
namespace NGS {

using namespace std;
/** \brief Default constructor, not PE and positive strand
 */
uTags::uTags():uGenericNGS()
{
}

/**< From uTokens */
uTags::uTags(const uToken & pToken)try:uGenericNGS(pToken){

if (pToken.isParamSet(token_param::CIGAR))
    setCigar(pToken.getParam(token_param::CIGAR));
 if (pToken.isParamSet(token_param::MAP_SCORE))
        setMapQual(std::stoi(pToken.getParam(token_param::MAP_SCORE)));
 if (pToken.isParamSet(token_param::PHRED_SCORE))
        setPhred(pToken.getParam(token_param::PHRED_SCORE));
 if (pToken.isParamSet(token_param::SEQUENCE))
        setSequence(pToken.getParam(token_param::SEQUENCE));
 if (pToken.isParamSet(token_param::FLAGS))
        setFlag(std::stoi(pToken.getParam(token_param::FLAGS)));
}
catch(ugene_exception_base &e)
{
        #ifdef DEBUG
        std::cerr << "Error in uGenericNGS(uToken)." <<std::endl;
        #endif
        e<<tag_error(*this);
        throw e;
}


/** \brief Copy constructor, with init list
 * \param otherItem: uGenericsNGS (or child class) object
 */
uTags::uTags(uGenericNGS otherItem):uGenericNGS(otherItem),name(nullptr),phredScore(nullptr),cigar(nullptr)
{
}

/** \brief Default constructor with init list, implicitly sets strand
 *
 * \param chr: name of the chromosome
 * \param start: beginning position of the tag
 * \param end: ending position of the tag
 */
uTags::uTags(std::string pchr, int start, int end, StrandDir pstrand):name(nullptr),phredScore(nullptr),cigar(nullptr)
{
   try {
     setStartEnd(start,end);
     setChr(pchr);
     setStrand(pstrand);
    }
    catch(elem_throw & e)
    {
        #ifdef DEBUG
               cerr << "Throwing in uTags constructor" <<endl;
        #endif

        string trace;
        if (std::string const * ste =boost::get_error_info<string_error>(e) )
                trace=*ste;

        e << string_error(trace+"Failling in uTags constructor, parameters are"+pchr+" "+std::to_string(start)+" "+std::to_string(end)+"\n");

        throw e;
    }
}

/** \brief Default destructor
 */
uTags::~uTags()
{

    if (name!=nullptr) {
        delete []name;
        name = nullptr;
    }

    if (phredScore!=nullptr) {
        delete []phredScore;
        phredScore = nullptr;
    }

    if (cigar!=nullptr) {
        delete []cigar;
        cigar = nullptr;
    }
}

//TODO TEST THIS
/** \brief Copy Constructor for UTags, necessary due to our char*
 * \param UTags& const copy_from: tag to copy
 */
uTags::uTags(const uTags& copy_from):uGenericNGS(copy_from),name(nullptr),phredScore(nullptr),cigar(nullptr)
{
    try {
    if (copy_from.name!=nullptr)
    {
        name = new char[(strlen(copy_from.name)+1)];
        strcpy(name,copy_from.name);
    }

    if (copy_from.cigar!=nullptr)
    {
        cigar = new char[(strlen(copy_from.cigar)+1)];
        strcpy(cigar,copy_from.cigar);
    }

    if (copy_from.phredScore!=nullptr)
    {
        phredScore = new char[(strlen(copy_from.phredScore)+1)];
        strcpy(phredScore,copy_from.phredScore);
    }
    strand=copy_from.strand;
    PELenght=copy_from.PELenght;
    sequence=copy_from.sequence;
    mapScore=copy_from.mapScore;
    Unmapped= copy_from.Unmapped;
    }
    catch(elem_throw & e){
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
    {throw e;};
}

/** \brief Overloaded assignement operator
 * \param uTags const& source
 * \return uTags& destination
 */
uTags& uTags::operator= (uTags const& assign_from)
{
    if (this == &assign_from) return *this;

    uGenericNGS::operator= (assign_from);

    if (name!=nullptr) {
        delete []name;
        name = nullptr;
    }

    if (phredScore!=nullptr) {
        delete []phredScore;
        phredScore = nullptr;
    }

    if (cigar!=nullptr) {
        delete []cigar;
        cigar = nullptr;
    }

    if (assign_from.name!=nullptr)
    {
        name = new char[(strlen(assign_from.name)+1)];
        strcpy(name,assign_from.name);
    }

    if (assign_from.cigar!=nullptr)
    {
        cigar = new char[(strlen(assign_from.cigar)+1)];
        strcpy(cigar,assign_from.cigar);
    }

    if (assign_from.phredScore!=nullptr)
    {
        phredScore = new char[(strlen(assign_from.phredScore)+1)];
        strcpy(phredScore,assign_from.phredScore);
    }
    strand=assign_from.strand;
    PELenght=assign_from.PELenght;
    sequence=assign_from.sequence;
    mapScore=assign_from.mapScore;
     Unmapped= assign_from.Unmapped;
    return *this;
}

/** \brief Associate a tag with it's PE mate. Currently unused
*
* \param pTag: Pointer to another Tag element to associate with
*
*/
void uTags::setMate(uTags* pTag)
{
    pMate=pTag;
}

// TODO: Should this be in parser code?
/** \brief Load a tag from a string formated in Sam format and enter in member data
 *
 * \param samString: String to load from
 * \param minimal : if true, do not load certain information for space contraints
 *
 */
void uTags::loadfromSamString(std::string samString, bool minimal=false)
{
    *this=factory::makeTagfromSamString(samString,minimal);
}

// TODO: Move this to output class
/** \brief Write our tag in bed format, using our Tag Name/ID as the misc 4th column
 *
 * \param out : Output stream to use.
 *
 */
void uTags::writeBedToOuput(std::ostream &out) const
{

    uGenericNGS::writeBedToOuput(out, false);
    if (getName().size()!=0)
        out << "\t" << getName();
    else
        out << "\t" << "unknown";

    out << "\t" << -1;

    /**< Strand is implicit, always write it */
    //if ((getStrand()=='+')||(getStrand()=='-'))
    char strand='+';
    if (isReverse())
        strand='-';
    out << "\t" << strand;

    out <<"\n";
}


// TODO: Move this to output class
/** \brief Write our tag in bed format, but only writing the + side of a paired end tag. Ignore others
 *
 * \param out : Output stream to use.
 *
 */
void uTags::writetoBedCompletePE( std::ostream &out)
{

    if ((this->isPE())&&(this->getStrand()==StrandDir::FORWARD))
    {
        out << getChr() << "\t" << this->getStart() << "\t" << ( this->getStart()+this->getPeLenght()) <<endl;;
    }

}

// TODO: Move this to output class
/** \brief Write our data in Sam format, with some minor details missing
 *
 * \param out : Output stream to use.
 *
 */
void uTags::writeSamToOutput(std::ostream &out) const
{

    int peLenght;
    string cigar, sequence, phred;

    cigar = this->getCigar();
    if (cigar.size()==0)
        cigar=utility::concatStringInt("M", getLenght(), false);
    sequence = this->getSequence();
    phred= getPhred();
    if (phred.size()==0)
        phred="*";
    /**< Since Picard tools poorly validates a sequence lenght, we must pad with = if we did not store sequences */
    if (sequence.size()==0)
    {
        for (int i=0; i<getLenght(); i++ )
            sequence.push_back('=');
    }

    if (isPE())
    {
        peLenght= getPeLenght();
        /**< We store lenght as positive, but in Sam format it can be negative */
        if (flag&0x10)
            peLenght = (-peLenght);
    }
    else
        peLenght=0;

    out <<  getName() << "\t" << getFlag() << "\t" << getChr() << "\t" << getStart() << "\t" <<getMapQual()
    << "\t" <<cigar << "\t" << '*' << "\t" << 0 << "\t" << peLenght << "\t" ;

    out << sequence  << "\t"<< phred;
    out <<endl;
}
/** \brief Function used for debugging
 */
void uTags::debugElem() const
{
    using namespace utility;
    stringTocerr("Outputting elemn data");
    stringTocerr("Chrom "+getChr());
    stringTocerr("Start "+std::to_string(getStart()));
    stringTocerr("End " +std::to_string(getEnd()));
    stringTocerr("PELenght " +std::to_string(getPeLenght()));
    stringTocerr("Flag " +std::to_string(getFlag()));
}

// TODO: Move this to output class
/** \brief Write our PE tags as a sam tag that is the lenght of the fragment associated with our PE
 *   Does nothing if the tag is not PE
 *
 * \param out : Output stream to use.
 *
 */
void uTags::writeCompletedPESamToOutput(std::ostream &out)
{
/**< Write function, do not trust for diverse use. */
try {
    int peLenght;
    peLenght= getPeLenght();

      if (flag&0x10)
            peLenght = (-peLenght);

        //TODO, introducting try catch
        if ((flag&0x10)||(peLenght<0)){
           #ifdef DEBUG
               cerr << " negative Strand??"<<endl;
            #endif

            debugElem();
          //  utility::pause_input();
        }
        else /**< Do not write cigar for now */
        if (peLenght!=0)
        {
            out <<  getName() << "\t" << 0 << "\t" << getChr() << "\t" << getStart() << "\t" <<getMapQual()
            << "\t" << peLenght<<'M' << "\t" << '*' << "\t" << 0 << "\t" << 0 << "\t" ;
            for (int i=0; i<peLenght; i++ )
                out << '=' ;
            out  << "\t"<< '*';

            out <<endl;
        }
}
catch(std::exception & e ){
        #ifdef DEBUG
            cerr << " Catching in writeCompletedPESamToOutput"<<endl;
        #endif
    throw e;
}

}

// TODO: Move this to output class
/** \brief Write a tag in Sam specification, but trimming a certain size both left and right of it.
 *
 * \param out : Our output stream
 * \param left: Number of bp to trim from left ( 3' direction ) of tag
 * \param right : Number of bp to trim from right ( 5' direction ) of tag.
 *
 */
bool uTags::writeTrimmedSamToOutput(std::ostream &out, int left, int right)
{

    // int ourflag;
    if ((left + right)>=this->getLenght())
    {
        return false;
    }
    // ourflag=ourflag&
    if (getStrand()==StrandDir::REVERSE)
    {
        trimSites(right, left);
    }
    else
    {
        trimSites(left,right);
    }    /**< Do not write cigar for now */
    out <<  getName() << "\t" << 0 << "\t" << getChr() << "\t" << getStart() << "\t" << getMapQual()
    << "\t" <<   getLenght()<<'M' << "\t" << '*' << "\t" << 0 << "\t" << 0;// << "\t" << '*' << "\t" << '*';
    out  << "\t";
    for (int i=0; i<getLenght(); i++ )
        out << '=' ;
    out  << "\t"<< '*';


    out <<endl;
    return true;
}

// TODO: Move to parser?
namespace factory
{
uTags makeTagfromSamString(std::string samString, bool minimal)
{
    stringstream Infostream;
    int size,ourFlag;
    Infostream.str(samString);
    string ourChr, ourStart, ourEnd, ourName;
    string phre, cig, seq;
    string temp,tempname;
    /**< Read name */
try {
    if (!minimal)
    {
        Infostream>> tempname;
        if (tempname.find("/")!=string::npos)
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
    string substr;
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

 //   ourEnd=(std::stoi(ourStart)+(size-1));
    uTags returnTag(ourChr,std::stoi(ourStart),std::stoi(ourStart)+(size-1) );
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
         returnTag.setStrand(StrandDir::REVERSE);// ='-';
     else
         returnTag.setStrand(StrandDir::FORWARD);// ='+'; */
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
	catch(elem_throw & e)
	{
	    string trace;
	     #ifdef DEBUG
                 cerr << "catching in factory on elem_throw" <<endl;
        #endif

	  if (std::string const * ste =boost::get_error_info<string_error>(e) )
			trace=*ste;

	   e << string_error(trace+"Failling in factory:makeTagfromSamString constructor,  on string \n"+samString);
	  throw e;
	}
	catch(...)
	{
                #ifdef DEBUG
                     cerr << "catching in factory on general throw" <<endl;
        #endif

	    elem_throw e;
	    e << string_error("we threw in makeTagfromSamString trying the next string \n"+samString);
	    throw e;
	}

}
} // End of namespace factory

/** \brief Empty constructor for uTagsChrom
 *
 */
uTagsChrom::uTagsChrom():uGenericNGSChrom<uTags>()
{

}

/** \brief Constructor thats sets the chrom used by the class
 *
 */

uTagsChrom::uTagsChrom(std::string ourChrom) : uGenericNGSChrom<uTags>(ourChrom)
{

}

// TODO: is it complete?
/** \brief Copy constructor ( well not really ). Not complete?
 *
 * \param copyCop : The object to instaciate  from.
 */
uTagsChrom::uTagsChrom(uGenericNGSChrom<uTags> copyCop)
{
    VecSites=copyCop.returnVecData();
    chr= copyCop.getChr();
}

/** \brief Output the chrom to bed format
 *
 * \param out std::ostream&  : Outputstream, needs to be opened and validated
 * \return void
 *
 */
void uTagsChrom::outputBedFormat(std::ostream& out) const
{
   // std::vector<uTags>::iterator iterVec;
    for (auto iterVec = VecSites.begin(); iterVec != VecSites.end(); ++iterVec)
    {
        iterVec->writeBedToOuput(out);
    }
}

/** \brief Write all PE tags are complete bed regions
 *
 * \param out : Our output stream
 *
 */
void uTagsChrom::writetoBedCompletePE(std::ostream& out)
{
    std::vector<uTags>::iterator iterVec;

    for (iterVec = VecSites.begin(); iterVec != VecSites.end(); ++iterVec)
    {
        iterVec->writetoBedCompletePE(out);
    }
}

/** \brief Write every tag in Sam format, trimming the size of every tag
 *
 * \param out : Our output stream
 * \param left: Number of bp to trim from left ( 3' direction ) of tag
 * \param right : Number of bp to trim from right ( 5' direction ) of tag.
 *
 */
void uTagsChrom::writeTrimmedSamToOutput(std::ostream &out, int left, int right)
{
    int errorcount=0;
    std::vector<uTags>::iterator iterVec;
   // if (VecSites.size())
   //     out << "@SQ	SN:" << getChr() << "\t"<<"LN:"<<getChromSize() <<endl;
    for (iterVec = VecSites.begin(); iterVec != VecSites.end(); ++iterVec)
    {
        if (iterVec->writeTrimmedSamToOutput(out, left, right)==false)
            errorcount++;
    }
    #ifdef DEBUG
    if (errorcount)
        cerr << "Skipped " << errorcount << " tags as trim would reduce them to 0. On Chrom " << this->getChr();
    #endif

}

/** \brief Write our PE tags as a sam file and complete
 *   Does nothing if the tag is not PE
 *
 *  \param out : Output stream to use.
 *
 */
void uTagsChrom::writeCompletedPESamToOutput(std::ostream &out)
{

    std::vector<uTags>::iterator iterVec;

/**< Do not output unmapped tags */

    for (iterVec = VecSites.begin(); iterVec != VecSites.end(); ++iterVec)
    {
        if ((iterVec->isPE())&&(iterVec->getStrand()==StrandDir::FORWARD)&&(iterVec->isMapped()))
            iterVec->writeCompletedPESamToOutput(out);
    }

}

/** \brief Write the @SQ line for the chromosome
 *
 * \param out std::ostream& Output stream to use
 * \return void
 *
 */
void uTagsChrom::writeSamHeaderLine(std::ostream &out) const
{
     if (VecSites.size())
        out << "@SQ	SN:" << getChr() << "\t"<<"LN:"<<getChromSize() <<endl;
}

//void writeSamToOutput(std::ostream &out);

/** \brief Output our data into Sam formati with legal headers
 *
 *  \param out : OfStream
 *
 */
void uTagsChrom::writeSamToOutput(std::ostream &out) const
{

    for (auto iterVec = VecSites.begin(); iterVec != VecSites.end(); ++iterVec)
    {
        iterVec->writeSamToOutput(out);
    }

}

/** \brief For PE data, find the mate of ever tag and associate there pointers. Structure subject to change
 *
 */
 //TODO remove this?
void uTagsChrom::matchPE()
{
    std::vector<uTags>::iterator iterVec;
    std::vector<uTags>::iterator innerVec;

    string current, next;

    for (iterVec = VecSites.begin() ; iterVec!= VecSites.end(); iterVec++)
    {

        if (iterVec->getMate()==NULL)
            for (innerVec = (iterVec +1) ; innerVec< VecSites.end(); ++innerVec)
            {
                current= iterVec->getName();
                next = innerVec->getName();
                if (iterVec->getName()==innerVec->getName() )
                {
                    /**< Assign pointer and complement */
                    iterVec->setMate(&(*innerVec));

                    innerVec->setMate(&(*iterVec));

                    break;
                }
            }
    }
}
//TODO Re-write this
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

    int pos=0;
    /**< Need to mess around with this later, make sure tags at the same position are not being messed with. */
    /**< We go one Kb earlier in our experiment, make sure we get all overlapping tags. */

    //TODO fix this
  //  pos = this->findPrecedingSite((start-1000), 0 , this->count()-1);

    /**< If no tag leftwise, we start at beginning */
    if (pos==-1)
        pos=0;

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

void uTagsExperiment::loadFromSamWithParser(std::string filepath)
{
    op_mode= ReadMode::DEFAULT;
    uParser ourParser(filepath, "SAM");
    auto chrList=ourParser.getHeaderParamVector(header_param::CHR);
    auto chrSizes=ourParser.getHeaderParamVector(header_param::CHR_SIZE);
try {
//        std::cerr << "name and size are :" <<chrList.at(i) <<" "<<chrSizes.at(i)<<std::endl;
    for (int i=0; i<(int)chrList.size();i++){
        this->setChromSize(chrList.at(i), std::stoi(chrSizes.at(i)));
    }

    while (ourParser.eof()==false){
        auto Token =ourParser.getNextEntry();
        this->addSite(uTags(Token));
    }
}
    catch(...){
        throw;
    }
}


/** \brief Loads an entire Sam file, SE or PE
 *
 * \param curStream : our input stream that should already be opened and validated
 * \param minimal : If true, do not load every data
 *
 */
void uTagsExperiment::loadFromSam(std::ifstream& curStream, bool minimal)
{
    string lineString;
    int count=0;
    ourStream = &curStream;
    op_mode= ReadMode::DEFAULT;
/**< Parse Header if header there is */
    parseSamHeader();
    try
    {
        while(!ourStream->eof())
        {
            getline(*ourStream, lineString);
            /**< Make sure this is not a header */
            if (ourStream->eof())
            {
                cerr<< "Finished loading" << count << " tags" <<endl;
                break;
            }
            {

                if (lineString.size()>5)
                    addSite(move(factory::makeTagfromSamString(lineString,minimal) ));
                else{
                #ifdef DEBUG
                    cerr <<"Skipping the following line as it does not satisfy minimal sam string size" <<endl;
                    cerr << lineString <<endl;
                #endif


                }
                count++;
            }
        }
    }
    catch(elem_throw & e)
     {
            #ifdef DEBUG
            cerr << "caught an exception in loadfromSam"<<endl;
                #endif
         string trace;

        if (std::string const * ste =boost::get_error_info<string_error>(e) )
                trace=*ste;

        e <<string_error(trace+"\n"+"falling from loadFromSam(stream, bool) while loading tag number"+std::to_string(count) );
            #ifdef DEBUG
                    cerr << "Throwing elem_throw" <<endl;
                #endif

        throw e;
     }
    catch(std::exception our_exception)
    {
          #ifdef DEBUG
                   cerr << "caught an exception in loadfromSam"<<endl;
                #endif
        throw;
    }
}

/** \brief Parses the header of a Sam file. Used for progressive loading
 *
 * \param curStream : Our input stream, should be opened and validated
 *
 */
void uTagsExperiment::loadSamHeader(std::ifstream& curStream)
{
    setFileStream(curStream);
    /**< Parse Header if header there is */
    parseSamHeader();
}

//TODO use parser
/** \brief Returns the tag parsed from the next line of our sam file. Abort if not in progressive mode
 *
 * \param minimal : If yes, we will not load everything from the line
 *
 */
uTags uTagsExperiment::nextSamLine(bool minimal)
{
    string lineString;
    uTags tempTag;

    if (op_mode==ReadMode::DEFAULT)
    {

        #ifdef DEBUG
            cerr << "Stream not in progressive mode" <<endl;
        #endif


        abort();
    }
    if (ourStream==NULL)
    {
         #ifdef DEBUG
             cerr << "Error, trying to load tag from empty stream"<<endl;
        #endif
        abort();
    }
    if (ourStream->peek()=='@')
    {
        #ifdef DEBUG
              cerr << "Error, header line found. Did you parse the headers?"<<endl;
        #endif
        abort();
    }
    getline(*ourStream, lineString);
    if (ourStream->eof()==false)
        tempTag.loadfromSamString(lineString,minimal);

    return tempTag;
}


/**< To be used when we want to read in gradual mode. Will read every header line and created the appropriate chromosome size */
/**< Objects */
/** \brief Parse the Sam header. Use for progressive loading, once our input stream is set.
 *
 */
void uTagsExperiment::parseSamHeader()
{
    string lineString;

    /*   if (op_mode==DEFAULT)
       {
           cerr << "Stream not in progressive mode, aborting from parseSamHeader()" <<endl;
           abort();
       }
    */
    if (ourStream==NULL)
    {

        #ifdef DEBUG
              cerr << "Error, trying to load tag from empty stream, aborting from parseSamHeader()"<<endl;
        #endif


        abort();
    }

    /**< If header */
    while (ourStream->peek()=='@')
    {
        getline(*ourStream, lineString);

        /**< Parse the line */
        stringstream Infostream;
        string temp,chrom, size;
        int chromsize;
        Infostream.str(lineString);
        /**< validate @SQ */
        Infostream >>temp;
        /**< Chrom size line? */

        if (temp.find("@SQ")!=string::npos)
        {
            string data;
            string chrom;
            while (!(Infostream.eof())){
                Infostream >> data;
                if (data.find("SN:")!=string::npos){
                    data.erase(0,3);
                    chrom=data;
                }
                if (data.find("AS:")!=string::npos){
                    /**< Skip */
                }
                if (data.find("LN:")!=string::npos){
                    data.erase(0,3);
                    chromsize=std::stoi(data);

                }

            }
             this->setChromSize(chrom, chromsize);
             //this->setChr
        }
    }
}

/** \brief Write our data in bed format
 *
 * \param out : Our output stream
 *
 */
void uTagsExperiment::writeToBed(ostream& out) const
{

    for (auto iterMap = ExpMap.begin(); iterMap != ExpMap.end(); ++iterMap)
    {
        iterMap->second.outputBedFormat(out);
    }
}

/** \brief Write our data in bed format, completing our PE tags and ignore lone mates
 *
 * \param out : Our output stream
 *
 */
void uTagsExperiment::writetoBedCompletePE(std::ostream& out)
{
    std::map<std::string,uTagsChrom>::iterator iterMap;

    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); ++iterMap)
    {
        iterMap->second.writetoBedCompletePE(out);
    }
}

/** \brief Write the entire dataset in completed PE Sam format
 *
 * \param out : Our output stream
 *
 */
void uTagsExperiment::writeCompletedPESamToOutput(std::ostream &out)
{
    for (auto iterMap = ExpMap.begin(); iterMap != ExpMap.end(); ++iterMap)
    {
        iterMap->second.writeSamHeaderLine(out);
    }
    //ChipPeakVectorMap::iterator iterMap;
    for (auto iterMap = ExpMap.begin(); iterMap != ExpMap.end(); ++iterMap)
    {
        iterMap->second.writeCompletedPESamToOutput(out);
    }
}

/** \brief Write our entire experiment in Sam format
 *
 * \param out std::ostream& : our Output stream
 * \return void
 */
void  uTagsExperiment::writeSamToOutput(std::ostream &out) const
{
   // std::map<std::string,uTagsChrom>::iterator iterMap;
    //ChipPeakVectorMap::iterator iterMap;

    for (auto iterMap = ExpMap.begin(); iterMap != ExpMap.end(); ++iterMap)
    {
        iterMap->second.writeSamHeaderLine(out);
    }

    for (auto iterMap = ExpMap.begin(); iterMap != ExpMap.end(); ++iterMap)
    {
        iterMap->second.writeSamToOutput(out);
    }
}

/** \brief Set the lenght limit of a chromosome
 *
 * \param chrom string : The unique ID of our chrom
 * \param size int : Max limit to set
 * \return void
 *
 */
void uTagsExperiment::setChromSize(string chrom, int size)
{
    uTagsChrom* ptempChrom;
    ptempChrom=&(ExpMap[chrom]);
    ptempChrom->setChr(chrom);
    ptempChrom->setChromSize(size);
}
//TODO set where this should go?
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

/** \brief Write our exp to sam output, trimming both sides.
 *
 * \param out std::ostream& : Our output stream
 * \param left int : Left trim
 * \param right int : Right trim
 * \return void
 *
 */
void uTagsExperiment::writeTrimmedSamToOutput(std::ostream &out, int left, int right)
{
    std::map<std::string,uTagsChrom>::iterator iterMap;

    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); ++iterMap)
    {
        iterMap->second.writeSamHeaderLine(out);
    }
    //ChipPeakVectorMap::iterator iterMap;
    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); ++iterMap)
    {
        iterMap->second.writeTrimmedSamToOutput(out, left, right);
    }
}
} // End of namespace NGS
