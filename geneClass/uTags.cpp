#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include "uTags.h"


using namespace std;

/** \brief Default constructor, not PE and positive strand
 */
uTags::uTags():uGenericNGS(),name(nullptr),phredScore(nullptr),cigar(nullptr)
{
    strand='+';

}

uTags::uTags(uGenericNGS otherItem):uGenericNGS(otherItem),name(nullptr),phredScore(nullptr),cigar(nullptr)
{
    strand='+';
}


/** \brief Default constructor with init list, implicitly sets strand
 *
 * \param
 * \param
 * \return
 *
 */
uTags::uTags(std::string chr, int start, int end):name(nullptr),phredScore(nullptr),cigar(nullptr)
{
    /**< We handle parent init, as we do not have the same sanity assumptions as the parent */
    /**< If start is after end, we assumed this is a negative tag */
     //TODO Remove
   try {
    if (start>end)
    {
        setStartEnd(end,start);
        setStrand('-');
    }
    else
    {
        setStartEnd(start,end);
        setStrand('+');
    }
    setChr(chr);

    }
    catch(elem_throw & e)
    {
        cerr << "Throwing in uTags constructor" <<endl;
        string trace;
        if (std::string const * ste =boost::get_error_info<string_error>(e) )
                trace=*ste;

        e << string_error(trace+"Failling in uTags constructor, parameters are"+chr+" "+utility::numberToString(start)+" "+utility::numberToString(end)+"\n");

        throw e;

    }

//    PE=false;
}

uTags::~uTags()    // default destructor
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

/** \brief Copy Constructor for UTags, necessary due to our char*
 *
 * \param UTags& const
 *
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
    catch(...){
         cerr << " Failling in copy constructor, uTags(uTags&) " <<endl;
    }
}

/** \brief Overloaded assignement operator
 *
 * \param
 * \param
 * \return
 *
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
    out << "\t" << getStrand();

    out <<"\n";
}

/** \brief Write our tag in bed format, but only writing the + side of a paired end tag. Ignore others
 *
 * \param out : Output stream to use.
 *
 */
void uTags::writetoBedCompletePE( std::ostream &out)
{

    if ((this->isPE())&&(this->getStrand()=='+'))
    {
        out << getChr() << "\t" << this->getStart() << "\t" << ( this->getStart()+this->getPeLenght()) <<endl;;
    }

}

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
    //Since Picard tools poorly validates a sequence lenght, we must pad with = if we did not store sequences
    if (sequence.size()==0)
    {
        for (int i=0; i<getLenght(); i++ )
            sequence.push_back('=');
    }

    if (isPE())
    {
        peLenght= getPeLenght();
        //We store lenght as positive, but in Sam format it can be negative
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

void uTags::debugElem() const
{
    using namespace utility;
    stringTocerr("Outputting elemn data");
    stringTocerr("Chrom "+getChr());
    stringTocerr("Start "+numberToString((int)getStart()));
    stringTocerr("End " +numberToString((int)getEnd()));
    stringTocerr("PELenght " +numberToString((int)getPeLenght()));
    stringTocerr("Flag " +numberToString((int)getFlag()));
}

/** \brief Write our PE tags as a sam tag that is the lenght of the fragment associated with our PE
 *   Does nothing if the tag is not PE
 *
 * \param out : Output stream to use.
 *
 */

//Write function, do not trust for diverse use.
void uTags::writeCompletedPESamToOutput(std::ostream &out)
{

    int peLenght;
    peLenght= getPeLenght();

      if (flag&0x10)
            peLenght = (-peLenght);

        //Todo, introducting try catch
        if ((flag&0x10)||(peLenght<0)){
           cerr << " negative Strand??"<<endl;
            debugElem();
          //  utility::pause_input();
        }
        else //Do not write cigar for now
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
/** \brief Write a tag in Sam specification, but trimming a certain size both left and right of it.
 *
 * \param out : Our output stream
 * \param left: Number of bp to trim from left ( 3' direction ) of tag
 * \param right : Number of bp to trim from right ( 5' direction ) of tag.
 *
 */

bool uTags::writeTrimmedSamToOutput(std::ostream &out, int left, int right)
{

    //  int ourflag;
    if ((left + right)>=this->getLenght())
    {
        return false;
    }
    // ourflag=ourflag&
    if (getStrand()=='-')
    {
        trimSites(right, left);
    }
    else
    {
        trimSites(left,right);
    }    //Do not write cigar for now
    out <<  getName() << "\t" << 0 << "\t" << getChr() << "\t" << getStart() << "\t" << getMapQual()
    << "\t" <<   getLenght()<<'M' << "\t" << '*' << "\t" << 0 << "\t" << 0;// << "\t" << '*' << "\t" << '*';
    out  << "\t";
    for (int i=0; i<getLenght(); i++ )
        out << '=' ;
    out  << "\t"<< '*';


    out <<endl;
    return true;
}

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
    //Read name
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

    //Read flag
    Infostream>>ourFlag;
    Infostream>>ourChr;

    //returnTag.setChr(ourChr);

    Infostream>>ourStart;

    int mapScore;
    Infostream>>mapScore;
    //Cigar
    Infostream>>cig;

    // Parse Cigar here
    //Find first letter , read number before, repeat
    //Easier with Regex,  but what the hell
    int curPos=0;
    string substr;
    size=0;
    for (unsigned int i=0; i< cig.size(); i++)
    {
        //If isAlpha then check previous numbers
        if (isalpha(cig.at(i)))
        {
            //If a count value
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

 //   ourEnd=(utility::stringToInt(ourStart)+(size-1));
    uTags returnTag(ourChr,utility::stringToInt(ourStart),utility::stringToInt(ourStart)+(size-1) );
    returnTag.setName(ourName);
    returnTag.setFlag(ourFlag);
    returnTag.setMapQual(mapScore);
    //name of next mate
    Infostream>>temp;
    //Pos of next mate
    Infostream>>temp;
    //Template lenght for PE readas
    int PELenght;
    Infostream>>PELenght;
    PELenght= abs(PELenght);
   // cerr << "PeLenght is" <<
    //Sequence
    Infostream>>seq;
    //Pred score of every position.
    Infostream>>phre;
    //Sam flag
    //StrongType this
     if (ourFlag&0x10)
         returnTag.setStrand('-');// ='-';
     else
         returnTag.setStrand('+');// ='+'; */
    if (ourFlag&0x4)
        returnTag.setMapped(false);
    else
        returnTag.setMapped(true);

    //PE validation
    //If PE and mate is aligned
    if ((ourFlag&0x1)&&(ourFlag&0x2))
    {
        returnTag.setPELenght(PELenght);
    }
    else
    {
        returnTag.setPELenght(0);
    }
//    pMate = NULL;
    //if we want to keep, we store, otherwise scope will erase our data
    if (!minimal)
    {
        returnTag.setPhred(phre) ;

        returnTag.setCigar(cig);
    }
    //Move semantics
    return returnTag;

}
catch(elem_throw & e)
{
    string trace;
    cerr << "catching in factory on elem_throw" <<endl;
  if (std::string const * ste =boost::get_error_info<string_error>(e) )
                trace=*ste;

   e << string_error(trace+"Failling in factory:makeTagfromSamString constructor,  on string \n"+samString);
  throw e;
}
catch(...)
{
      cerr << "catching in factory on general throw" <<endl;
    elem_throw e;
    e << string_error("we threw in makeTagfromSamString trying the next string \n"+samString);
    throw e;
}

}
}

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

/** \brief Copy constructor ( well not really ). Not complete?
 *
 * \param copyCop : The object to instaciate  from.
 * \param
 * \return
 *
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

    if (errorcount)
        cerr << "Skipped " << errorcount << " tags as trim would reduce them to 0. On Chrom " << this->getChr();
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
        if ((iterVec->isPE())&&(iterVec->getStrand()=='+')&&(iterVec->isMapped()))
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
 *
 */

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
                    //Assign pointer and complement
                    iterVec->setMate(&(*innerVec));

                    innerVec->setMate(&(*iterVec));

                    break;
                }
            }
    }
}
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

    int pos;
    //Need to mess around with this later, make sure tags at the same position are not being messed with.
    //We go one Kb earlier in our experiment, make sure we get all overlapping tags.
    pos = this->findPrecedingSite((start-1000), 0 , this->count()-1);

    //If no tag leftwise, we start at beginning
    if (pos==-1)
        pos=0;

    iterVec=VecSites.begin();

    for (iterVec=(iterVec+pos) ; iterVec != VecSites.end(); ++iterVec)
    {
        //If our start is passed our region  we stop. time saver.
        if (iterVec->getStart()> end)
            break;


        if (overlap==false)
        {

            if (utility::isRegionAInsideRegionB(iterVec->getStart(), iterVec->getEnd(), start, end))
            {
                //Increase signal for every overlapping position of the tag.
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
    op_mode= DEFAULT;
//Parse Header if header there is
    parseSamHeader();
    try
    {
        while(!ourStream->eof())
        {
            getline(*ourStream, lineString);
            //Make sure this is not a header
            if (ourStream->eof())
            {
                cerr<< "Finished loading" << count << " tags" <<endl;
                break;
            }
            {
                //  uTags tempTag;
                //  tempTag.loadfromSamString(lineString, minimal);
                //tempPeak.loadfromBedString(peakString);

                if (lineString.size()>5)
                    addSite(move(factory::makeTagfromSamString(lineString,minimal) ));
                else{
                    cerr <<"Skipping the following line as it does not satisfy minimal sam string size" <<endl;
                    cerr << lineString <<endl;
                }
                count++;
                if (count%2000000==0)
                    cerr << count <<endl;
            }
        }
    }
    catch(elem_throw & e)
     {
         cerr << "caught an exception in loadfromSam"<<endl;
         string trace;

        if (std::string const * ste =boost::get_error_info<string_error>(e) )
                trace=*ste;

        e <<string_error(trace+"\n"+"falling from loadFromSam(stream, bool) while loading tag number"+utility::numberToString(count) );

        cerr << "Throwing elem_throw" <<endl;
        throw e;
     }
    catch(std::exception our_exception)
    {
        cerr << "caught an exception in loadfromSam"<<endl;
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


/** \brief Returns the tag parsed from the next line of our sam file. Abort if not in progressive mode
 *
 * \param minimal : If yes, we will not load everything from the line
 *
 */
uTags uTagsExperiment::nextSamLine(bool minimal)
{
    string lineString;
    uTags tempTag;

    if (op_mode==DEFAULT)
    {
        cerr << "Stream not in progressive mode" <<endl;
        abort();
    }
    if (ourStream==NULL)
    {
        cerr << "Error, trying to load tag from empty stream"<<endl;
        abort();
    }
    if (ourStream->peek()=='@')
    {
        cerr << "Error, header line found. Did you parse the headers?"<<endl;
        abort();
    }
    getline(*ourStream, lineString);
    if (ourStream->eof()==false)
        tempTag.loadfromSamString(lineString,minimal);

    return tempTag;
}


//To be used when we want to read in gradual mode. Will read every header line and created the appropriate chromosome size
//Objects

/** \brief Parse the Sam header. Use for progressive loading, once our input stream is set.
 *
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
        cerr << "Error, trying to load tag from empty stream, aborting from parseSamHeader()"<<endl;
        abort();
    }

    //If header
    while (ourStream->peek()=='@')
    {
        getline(*ourStream, lineString);

        //Parse the line
        stringstream Infostream;
        string temp,chrom, size;
        int chromsize;
        Infostream.str(lineString);
        //validate @SQ
        Infostream >>temp;
        //Chrom size line?

        if (temp.find("@SQ")!=string::npos)
        {
            string data;
            string chrom;
            while (!(Infostream.eof())){
                Infostream >> data;
                if (data.find("SN:")!=string::npos){
                    data.erase(0,3);
                    chrom=data;
                  //  cerr << " Chrom is " << chrom << endl;
                }
                if (data.find("AS:")!=string::npos){
                   // data.erase(0,3);
                  //  chrom=data;

                }
                if (data.find("LN:")!=string::npos){
                    data.erase(0,3);
                    chromsize=utility::stringToInt(data);
                   // cerr << " ChromSize is " << chromsize << endl;
                }


               // cerr << data <<endl;
                               /* //Get Chrom
                Infostream >>chrom;
                //Erase SN:
                chrom.erase(0,3);

    //TODO repair , not parsing correctly

                Infostream>>size;
                //Erase LN:
                size.erase(0,3);
                chromsize=utility::stringToInt(size);
                 */
            }
             this->setChromSize(chrom, chromsize);
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
 * \param
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
    ptempChrom->setChromSize(size);
}

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


