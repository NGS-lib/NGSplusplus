#ifndef UFORMATEXPERIMENT_H_INCLUDED
#define UFORMATEXPERIMENT_H_INCLUDED
#include <fstream>
#include "IO/Parser/uParser.h"
#include <functional>
#include "uFormatChrom.h"

namespace NGS
{

// TODO: ComputeOnAllSites

enum class ReadMode
{
    DEFAULT, GRADUAL
};

//NOTE,Discuss use of a function that returns all chrom names.
/**<
Explanation on _SELF_
This is a simple implementation of the curiously recurring template (CRT)
Note that applying it to our entire class does contribute to code bloat
as functions that do not need to know their type, will still generate a seperate instantiation
for each data type.

Potentially, we could derive a parent class where we place the version that do not require _SELF_

However, while this would reduce code size, it would complexify an already somewhat complex hiearchy.
_BASE_ is our Tags, _CHROM_ our Chrom structure.

 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
class uGenericNGSExperiment
{
    /**< Can be used to make sure EXP, chrom both use appropriate base */
    static_assert(
        std::is_convertible<_CHROM_, uGenericNGSChrom<_CHROM_,_BASE_>>::value,
        "Both types do not use the same underlying data structures"
    );

    /**< Iterator typedef */
    typedef std::map<std::string,_CHROM_>      NGSExpMap;
    typedef typename NGSExpMap::iterator       NGSExpIter;
    typedef typename NGSExpMap::const_iterator NGSExpConstIter;
    typedef typename NGSExpMap::value_type     NGSExpPair;

    typedef std::vector<_BASE_> VecGenericNGS;
    typedef typename std::vector<_BASE_>::iterator VecGenIter;
    typedef typename std::vector<_BASE_>::const_iterator VecGenConstIter;

private:

    /**< Comparison functors */
    static bool comparePosStart(const _BASE_ &item, const int & value)
    {
        return item.getStart() < value;
    }
    static bool compareStart(const _BASE_ &item1, const _BASE_ &item2)
    {
        return item1.getStart() < item2.getStart();
    }
    static bool compareLenght(const _BASE_ &item1, const _BASE_ &item2)
    {
        return item1.getLenght() < item2.getLenght();
    }
    static bool comparePos(const _BASE_ &item1, const _BASE_ &item2)
    {
        if ((item1.getStart() != item2.getStart()))
            return item1.getStart() < item2.getStart();
        else
            return item1.getEnd() < item2.getEnd();
    }


protected:

    std::map<std::string,_CHROM_>  ExpMap= {};

    std::function<float(const _BASE_*)> sortGetStart=nullptr;
    std::function<float(const _BASE_*)> sortGetEnd=nullptr ;
    std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> m_comptFunc=compareStart;

    void removeSite(const std::string & pChr,const long int position);
    void removeSite(const std::string & pChr,const long int pStart,const long int pEnd);


public:
    /** \brief Empty constructor. Does nothing.
      */
    virtual ~uGenericNGSExperiment() {};
    uGenericNGSExperiment& operator=(const uGenericNGSExperiment& copFrom)=default;
    uGenericNGSExperiment(const uGenericNGSExperiment&)=default;

    virtual void addData(const _BASE_ &);
    virtual void addData(const _CHROM_ &);
    virtual void addData(const uToken &);

    void replaceChr(const _CHROM_ &);
    void removeChr(const std::string &);
    bool isChrom(const std::string &)const;
    long long count() const;
    void sortSites();
    template<typename Compare>
    void sortSites(Compare comp,std::function<float(const _BASE_*)> getStart_funct=nullptr,std::function<float(const _BASE_*)> getEnd_funct=nullptr);
    bool isSorted()const;
    typename std::vector<_BASE_>::const_iterator findPrecedingSite(std::string chr, int position)const;
    typename std::vector<_BASE_>::const_iterator findNextSite(std::string chr, int position)const;

    void removeSite(VecGenConstIter pItrPos);
    void removeSite(VecGenConstIter pItrStart,VecGenConstIter pItrEnd);

    virtual void loadWithParser(std::ifstream&, std::string,long long=0);
    virtual void loadWithParser(std::string, std::string,long long=0);
    virtual void loadWithParser(uParser&, long long=0);

    template<class UnaryPredicate>
    void loadWithParser_if(uParser& pParser,UnaryPredicate predicate, long long pBlockCount=0);

    template<class UnaryFunction>
    void loadWithParserAndRun(std::ifstream& pStream, std::string pType, UnaryFunction funct , int pBlockSize=1);
    template<class UnaryFunction>
    void loadWithParserAndRun(std::string filepath, std::string pType, UnaryFunction f, int pBlockSize=1);
    template<class UnaryFunction>
    void loadWithParserAndRun(uParser& pParser, UnaryFunction funct , int pBlockSize=1);

    void inferChrSize();
    void writeWithWriter(uWriter& pWriter) const;

    auto begin()->decltype(ExpMap.begin())
    {
        return ExpMap.begin();
    };
    auto end()->decltype(ExpMap.end())
    {
        return ExpMap.end();
    };

    auto begin()const->decltype(ExpMap.cbegin())
    {
        return ExpMap.cbegin();
    };
    auto end()const->decltype(ExpMap.cend())
    {
        return ExpMap.cend();
    };



    bool F(const std::string & pChrom) const;

    _CHROM_ getChrom(const std::string & chrom) const;
    const _CHROM_* getpChrom(const std::string & chrom) const;
    _CHROM_* getpChrom(const std::string & chrom);
    _BASE_ getSite(const std::string & pChr, long int pPosition)const;
    _BASE_ getSite(typename std::vector<_BASE_>::const_iterator posItr)const;

    _SELF_ getOverlapping(_SELF_ &compareExp, OverlapType type=OverlapType::OVERLAP_PARTIAL);
    _SELF_ getOverlapping(_CHROM_ &compareChrom, OverlapType type=OverlapType::OVERLAP_PARTIAL);
    _SELF_ getOverlapping(std::string chr, int start, int end, OverlapType type=OverlapType::OVERLAP_PARTIAL);

    // TODO: getSubset overload that check every chromosome
    _CHROM_ getSubset(const std::string & pChr, const double pStart, const double pEnd, OverlapType options=OverlapType::OVERLAP_PARTIAL);
    _CHROM_ removeSubset(const std::string & pChr,const double pStart, const double pEnd, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);

    // TODO: getDistinct overload that check every chromosome
    // TODO: getDistinct overload with range vector
    _SELF_ getDistinct(const std::string & pChr, const double pStart, const double pEnd, OverlapType type=OverlapType::OVERLAP_PARTIAL);
    _SELF_ removeDistinct(const std::string & pChr,const double pStart, const double pEnd, OverlapType options=OverlapType::OVERLAP_PARTIAL);


    long int getSubsetCount(const std::string & pChr, const double pStart, const double pEnd, const OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
    long int getSubsetCount(const _BASE_ & subsetReg, const OverlapType overlap=OverlapType::OVERLAP_PARTIAL);


    int getChrCount(){return ExpMap.size();};

    void setChrSize(std::string chr, int chrSize);
    int getChrSize(std::string chr);
    void divideItemsIntoBinofSize(int N, SplitType type=SplitType::STRICT);
    void divideItemsIntoNBins(int N, SplitType type=SplitType::STRICT);

    template<class UnaryPredicate>
    void removeSpecificSites(UnaryPredicate pred);
    /**< Wrappers around the STL algorithms */
    template<class BinaryOperation, class InitialValue>
    InitialValue accumulateChromsInfo(BinaryOperation binary_op, InitialValue init) const;
    template<class UnaryOperation>
    auto computeOnAllChroms(UnaryOperation unary_op) const -> std::map<std::string, decltype(unary_op(_CHROM_()))>;

    template<class UnaryFunction>
    UnaryFunction applyOnAllChroms(UnaryFunction f);
    template<class UnaryFunction>
    UnaryFunction applyOnAllChroms(const UnaryFunction f)const;
    template<class UnaryFunction>
    UnaryFunction applyOnSites(UnaryFunction f);
    template<class UnaryFunction>
    UnaryFunction applyOnSites(const UnaryFunction f)const;

    template <class UnaryPredicate>
    typename std::iterator_traits<NGSExpIter>::difference_type
    countChromsWithProperty(UnaryPredicate pred) const;
    template<class Compare>
    NGSExpConstIter maxChrom(Compare comp) const;
    template<class Compare>
    NGSExpConstIter minChrom(Compare comp) const;
    template<class Compare>
    std::pair<NGSExpConstIter, NGSExpConstIter> minAndMaxChroms(Compare comp) const;
//    template<class _SELF_, typename _CHROM_, typename _BASE_>
    /**< End STL wrappers */


      /** \brief Get the chromosomes for which a certain predicate is true
      *
      *  Returns a subset of the chromosome structures that evaluated true to the given predicate.
      *  The passed predicated must return true or false.
      *
      *
      * \param p UnaryPredicate : Unary predicate to evaluate on all chromosomes
      * \return A collection containing all the chromosomes for which the predicate is true
      */
    //NOTE, should this return an experiment?
    template<class UnaryPredicate>
    auto getSpecificChroms(UnaryPredicate pred) const->decltype(ExpMap)
    {
        decltype(ExpMap) copyColl;
        copy_if(std::begin(ExpMap), std::end(ExpMap), std::inserter(copyColl, std::begin(copyColl)), [&pred]( const typename decltype(ExpMap)::value_type& element)
        {
            return pred(element.second);
        });
        return copyColl;
    }

    uGenericNGSExperiment() {};

};


/** \brief Add a site to the experiment.
 *
 * \param const _BASE_ & newSite: the site to add to the experiment.
 * \sa addData
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::addData(const _BASE_ & newSite)
{
    try
    {
        _CHROM_* ptempChrom;
        ptempChrom=&(ExpMap[newSite.getChr()]);
        ptempChrom->setChr(newSite.getChr());
        ptempChrom->addData(newSite);
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Catching and re-throwing in uFormatExp::addData()" <<std::endl;
#endif
        throw;
    }
}
template<class _SELF_, typename _CHROM_, typename _BASE_>
/** \brief Transform a token and add the necessary information. Typically this creates a single site
 *
 * \param pToken const uToken&
 * \sa addData
 * \return void uGenericNGSExperiment<_SELF_,_CHROM_,
 *
 */
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::addData(const uToken & pToken)
{
    try
    {
        std::string chr = pToken.getParam(token_param::CHR);
        _CHROM_* ptempChrom;
        ptempChrom=&(ExpMap[chr]);
        ptempChrom->setChr(chr);
        ptempChrom->addData(pToken);
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Catching and re-throwing in uFormatExp::addData()" <<std::endl;
#endif
        throw;
    }


}


/** \brief Remove a specific number from the specific scaffold. It is recommende to work through the iterator overloads
 *
 * \param chr std::string : Key to map element (chrom)
 * \param position int : position to remove
 * \return void uGenericNGSExperiment<_CHROM_,
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::removeSite(const std::string & pChr, const long int pPosition)
{
    _CHROM_* tempChrom;
    tempChrom=&(ExpMap[pChr]);
    try
    {
        if (pPosition>=tempChrom->count())
        {
            throw param_throw()<<string_error("Crashing in removeSite() in uGenericNGSEXperiment, required index higher then number of sites");
        }
        tempChrom->removeSite(pPosition);

    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Catching and re-throwing in removeSite() in uGenericNGSEXperiment,";
#endif
        throw;
    }
}

/** \brief Remove a specific range of elements from the specific scaffold. It is recommende to work through the iterator overloads
 *
 * \param chr std::string : Key to map element (chrom)
 * \param position int : position to remove
 * \return void uGenericNGSExperiment<_CHROM_,
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::removeSite(const std::string &  pChr, const long int pStart, const long int pEnd)
{
    try
    {
        _CHROM_* tempChrom;
        if (this->isChrom(pChr)==false)
            throw param_throw()<<string_error("Throwing in removeSite(long int,long int). Required scaffold that is not set");

        tempChrom=&(ExpMap[pChr]);
        if ((pStart>pEnd)||(!pStart)||(pEnd))
            throw param_throw()<<string_error("Throwing in removeSite(long int,long int). Either index under 0, or end position smaller then start");

        if (pEnd>=tempChrom->count())
        {
            throw param_throw()<<string_error("Crashing in removeSite() in uGenericNGSEXperiment, required index higher then number of sites");
        }
        tempChrom->removeSite(pStart,pEnd);

    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Crashing in removeSite() in uGenericNGSEXperiment,";
#endif
        throw;
    }
}


/** \brief Remove the element pointed at by iterator. This function is a relatively straightforward wrapping of
 *         std::erase and follows the same rules and canveats
 * \param chr std::string : Key to map element (chrom)
 * \param position int : position to remove
 * \return void uGenericNGSExperiment<_CHROM_,
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::removeSite(VecGenConstIter pItrPos)
{
    try
    {
        _CHROM_* tempChrom;
        if (this->isChrom(pItrPos->getChr())==false)
            throw param_throw()<<string_error("Throwing in removeSite(VecGenConstIter). Required scaffold that is not set");

        tempChrom=&(ExpMap[pItrPos->getChr()]);
        tempChrom->removeSite(pItrPos);
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Crashing in removeSite() in uGenericNGSEXperiment,";
#endif
        throw;
    }
}


/** \brief Remove the range of elements pointed at by the iterator. This function is a relatively straightforward wrapping of
 *         std::erase and follows the same rules and canveats. Reminder, std:: erase will not erase the element pointed by the end-of-range iterator.
 *
 * \param chr std::string : Key to map element (chrom)
 * \param start VecGenConstIter : Iterator pointing to the start of the range to erase
 * \param end VecGenConstIter : Iterator pointing to the end of the range to erase
 * \sa removeSite
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::removeSite(VecGenConstIter pItrStart,VecGenConstIter pItrEnd)
{
    try
    {
    	if (pItrStart->getChr() != pItrEnd->getChr()) {
		throw param_throw() << string_error("Throwing in removeSite(VecGenConstIter, VecGenConstIter). Iterators on different chromosomes.");
	}

        _CHROM_* tempChrom;
        if (this->isChrom(pItrStart->getChr())==false) {
            throw param_throw()<<string_error("Throwing in removeSite(VecGenConstIter, VecGenConstIter). Required scaffold that is not set");
	}

        tempChrom=&(ExpMap[pItrStart->getChr()]);
        tempChrom->removeSite(pItrStart,pItrEnd);
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Crashing in removeSite() in uGenericNGSEXperiment,";
#endif
        throw;
    }
}



/** \brief Infer and set the of every chrom base on their region at the highest position.
* \return void
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::inferChrSize()
{
    applyOnAllChroms(std::mem_fun_ref(static_cast<void (_CHROM_::*)()>(&_CHROM_::inferChrSize)));
    //applyOnAllChroms(std::bind(static_cast<void (_CHROM_::*)()>(&_CHROM_::inferChrSize)));
}


/** \brief load basic data from a Parser and load necessary data by passing to object constructor
 *
 * \param stream std::ifstream& file to load from
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::loadWithParser(std::ifstream& pStream, std::string pType, long long pBlockCount)
{
    try
    {
        std::istream& refStream = pStream;
        uParser Curparser(&refStream, pType,pBlockCount);
        loadWithParser(Curparser);
    }
    catch (...)
    {
        throw ;
    }
}

/** \brief Load a file
 * \param std::string filepath: the path to the file to load
 * \param std::string pType: the file type (i.e.: BED, SAM, BEDGRAPH, WIG, etc...)
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::loadWithParser(std::string filepath, std::string pType, long long pBlockCount)
{

    try
    {
        uParser ourParser(filepath, pType,pBlockCount);
        loadWithParser(ourParser);
    }
    catch (...)
    {
        throw;
    }
}

/** \brief Load a file
 * \param std::string filepath: the path to the file to load
 * \param std::string pType: the file type (i.e.: BED, SAM, BEDGRAPH, WIG, etc...)
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::loadWithParser(uParser& pParser, long long pBlockCount)
{
    try
    {
        if (pBlockCount>0)
        {
            int counter=0;
            while ((pParser.eof()==false)&&(counter!=pBlockCount))
            {
                this->addData(pParser.getNextEntry());
                counter++;
            }
        }
        else
        {
            {
                while (pParser.eof()==false)
                    this->addData((pParser.getNextEntry()));
            }
        }
    }
    catch (...)
    {
        throw;
    }
}

/** \brief Load up to pBlockCount data (0 = load entire file), keeping only data that fits condition
 * \param std::string filepath: the path to the file to load
 * \param std::string pType: the file type (i.e.: BED, SAM, BEDGRAPH, WIG, etc...)
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryPredicate>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::loadWithParser_if(uParser& pParser,UnaryPredicate predicate, long long pBlockCount)
{
    try
    {
        if (pBlockCount)
        {
            int counter=0;
            while ((pParser.eof()==false)&&(counter!=pBlockCount))
            {
                uToken loadedToken=pParser.getNextEntry();
                if (predicate(_BASE_(loadedToken)))
                {
                    this->addData(loadedToken);
                    counter++;
                }
            }
        }
        else
        {
            {
                while (pParser.eof()==false)
                    this->addData((pParser.getNextEntry()));
            }
        }
    }
    catch (...)
    {
        throw;
    }
}


/** \brief Write our data using the provided Writer
 *
 * \param out std::ofstream& stream to write to
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::writeWithWriter(uWriter& pWriter) const
{
    try {
    auto writeFunct= std::bind(&_CHROM_::writeWithWriter,std::placeholders::_1, std::ref(pWriter));
    applyOnAllChroms(writeFunct);
    }
    catch(...){throw;}
}


/** \brief Check if a chrom collection associated with the passed ID exist
 * \param const std::string & chrom: the name of the collection.
 * \return bool : True if the chrom collection exist
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
bool uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::isChrom(const std::string & pChrom) const
{
    try {
    return (ExpMap.count(pChrom));
    }
    catch(...){throw;}
}

/** \brief Returns the requested chrom object
 * \param const std::string & chrom: the name of the chrom.
 * \exception ugene_operation_throw: When the name of the chrom does not exists.
 * \return _CHROM_: a copy(?) of the chrom object
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_CHROM_ uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getChrom(const std::string & chrom) const
{   try {
    if (ExpMap.count(chrom)==0)
    {
        throw param_throw()<<string_error("Requested non-existent Chrom from Exp in getChrom(), value : " +chrom);
    }
    return ExpMap.find(chrom)->second.getCopy();
    }
    catch(...)
    {
        throw;
    }
}


/** \brief Returns a const pointer to the requested chrom object, if it exists.
 * \param const std::string & chrom: the name of the chrom.
 * \exception ugene_operation_throw: When the name of the chrom does not exists.
 * \return const _CHROM_*: a const pointer to the chrom object
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
const _CHROM_* uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getpChrom(const std::string & chrom) const
{
    try {
    if (ExpMap.count(chrom)==0)
    {
        throw param_throw()<<string_error("Required pointer to non-existent Chrom from Exp in getpChrom(), value : " +chrom);
    }
    const auto refer=&(ExpMap.find(chrom)->second);
    return (refer);
    }
    catch(...){throw;}
}

/** \brief Returns a pointer to the requested chrom object, if it exists.
 * \param const std::string & chrom: the name of the chrom.
 * \exception ugene_operation_throw: When the name of the chrom does not exists.
 * \return _CHROM_*: a pointer to the chrom object
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_CHROM_* uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getpChrom(const std::string & chrom)
{
    try {
    if (ExpMap.count(chrom)==0)
    {
        throw param_throw()<<string_error("Required pointer to non-existent Chrom from Exp in getpChrom(), value : " +chrom);
    }
    return &(ExpMap[chrom]);
    }
    catch(...){throw;}
}

/** \brief Set the size of a chrom object
 *
 *  Take note, the scaffold size is no guarantee. Various tools may map
 *  elements over the end of a reference/chr. As such, scaffold size is provided as is.
 *
 * \param std::string chr: the chrom from which to set the size.
 * \param int chrSize: the size of the chrom.
 * \return void
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::setChrSize(std::string chr, int chrSize)
{
    try
    {
        getpChrom(chr)->setChromSize(chrSize);
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Get the size of a chrom object
 * \param std::string chr: the chrom from which to set the size.
 * \return int: the size of the chromosome.
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
int uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getChrSize(std::string chr)
{
    try
    {
        return getpChrom(chr)->getChromSize();
    }
    catch(...)
    {
        throw;
    }
}

//Return the number of elements in our experiment
/** \brief Return our element count
 *
 * \return int : Number of elements
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
long long uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::count() const
{
    return accumulateChromsInfo([](long long partialSum, _CHROM_ chrom)
    {
        return partialSum + chrom.count();
    }, 0LL);
}

/** \brief Return overlap of a specific position for a specific map
 *
 * \param chr std::string : Map Key ( chrom)
 * \param start int : Start position, must be positive
 * \param end int : End position, must be >= start
 * \param overlap int : Type of overlap
 * \exception ugene_operation_throw: when the chr value is not valid
 * \return int: the count.
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
long int uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getSubsetCount(const std::string & pChr, const double pStart, const double pEnd, OverlapType overlap)
{
    long int count = 0;
    try
    {
        count = getpChrom(pChr)->getSubsetCount(pStart, pEnd, overlap);
    }
    catch (...)
    {
        throw;
    }
    return count;
}

/** \brief Return overlap of a specific position
 *
 * \param const _BASE_ & subsetReg: a tag that represent a subset of an experiment.
 * \param overlap int : Type of overlap
 * \return int: the count.
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
long int uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getSubsetCount(const _BASE_ & subsetReg, const OverlapType overlap)
{
    try
    {
        long int count=0;
        count = getSubsetCount(subsetReg.getChr(),
                               subsetReg.getStart(),
                               subsetReg.getEnd());
        return count;
    }
    catch(...)
    {
        throw;
    }
}


/** \brief Sort Every site of every chrom based on location
 *
 * \return void uGenericNGSExperiment<_CHROM_,
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::sortSites()
{
    sortSites(uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::compareStart,&_BASE_::getStart,&_BASE_::getEnd);
}

/** \brief Sort the sites vector by applying a certain comparison
    *      See documention on Chrom version for indepth comments
    *
    * \param comp Compare : Binary comparison operation to perform on the sites collection
    * \return void
    */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<typename Compare>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::sortSites(Compare comp,std::function<float(const _BASE_*)> getStart_funct,std::function<float(const _BASE_*)> getEnd_funct)
{
    try
    {
        sortGetStart=getStart_funct;
        if (getEnd_funct==nullptr)
            sortGetEnd=sortGetStart;
        else
            sortGetEnd= getEnd_funct;

        m_comptFunc=comp;
        /**< AS there are two SorSites functions, we must specify this rather clunky signature so it knows what overload to use */
        auto sortfunct=std::bind( (void(_CHROM_::*)(Compare,std::function<float(const _BASE_*)>,std::function<float(const _BASE_*)>))
                                  &_CHROM_::sortSites,std::placeholders::_1, comp,getStart_funct,getEnd_funct);
        applyOnAllChroms(sortfunct );
        // return std::sort(std::begin(VecSites), std::end(VecSites), comp);
    }
    catch(std::exception &e)
    {
        throw e;
    }
}


/** \brief Returns false if at least one chrom is unsorted
 *
 * \return bool true if the experiment is sorted
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
bool uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::isSorted()const
{
    try {
    bool sorted=true;
    applyOnAllChroms([&](const _CHROM_& chrom)
    {
        if (chrom.isSorted(m_comptFunc)==false)
            sorted=false;
    });

    return sorted;
    }
    catch(...){throw;}
}

/** \brief Return an interator pointing to the element of the chr before the specified value
 *   Note that this is based on the current sort type so may not refer to genomic position.
 *   Requires the data to be sorted first
 * \param std::string chr: chrom to search
 * \param int position: value to evaluate from (based on the sort type: i.e.: if sorted by score, it will use score as position)
 * \exception param_throw: When the chr value is not present in current experiment.
 * \return Iterator pointing to value.
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
typename std::vector<_BASE_>::const_iterator uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::findPrecedingSite(std::string chr, int position)const
{

    if (!ExpMap.count(chr))
    {
        throw param_throw()<<string_error("Failling in uGenericNGSExperiment::findPrecedingSite, value "+chr+" does not exist.\n");
    }
    try
    {

        auto tempChrom = getpChrom(chr);
        return tempChrom->findPrecedingSite(position);
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Return an interator pointing to the element of the chr after the specified value
 *   Note that this is based on the current sort type so may not refer to genomic position.
 *   Requires the data to be sorted first
 * \param std::string chr: chrom to search
 * \param int position: value to evaluate from (based on the sort type: i.e.: if sorted by score, it will use score as position)
 * \exception param_throw: When the chr value is not present in current experiment.
 * \return Iterator pointing to value
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
typename std::vector<_BASE_>::const_iterator uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::findNextSite(std::string chr, int position)const
{
    try
    {
        if (!ExpMap.count(chr))
        {
            throw param_throw()<<string_error("Failling in uGenericNGSExperiment::findNextSite, value "+chr+" does not exist.\n");
        }
        auto tempChrom = getpChrom(chr);
        return tempChrom->findNextSite(position);
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Get a specific site from a specific chrom. Overloaded to work with position, typically got from findPrecedingor findNext
 *
 * \param chr std::string
 * \param position int
 * \return _BASE_
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_BASE_ uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getSite(const std::string & pChr, const long int pPosition) const
{
    try
    {
        auto tempChrom = getpChrom(pChr);
        return tempChrom->getSite(pPosition);
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Get a specific site from a specific chrom. Overloaded to work with an interator, typically got from findPrecedingor findNext
 *
 * \param chr std::string
 * \param position int
 * \return _BASE_
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_BASE_ uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getSite(typename std::vector<_BASE_>::const_iterator posItr)const
{
    return *posItr;
}

/** \brief Return a Chrom containing only the sites that overlap the given chr. REQUIRES COLLECTION TO BE SORTED
 *
 * \param chr std::string : Name of scaffold to subset on
 * \param start int : Start position
 * \param end int : End position
 * \param options OverlapType : OverlapType if needed
 * \return _CHROM_ : Chrom containing the overlapping elements
 *
 */
 // TODO: should return a sorted chrom
template<class _SELF_, typename _CHROM_, typename _BASE_>
_CHROM_ uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getSubset(const std::string & pChr, const double pStart, const double pEnd, OverlapType options)
{
    try {
    if (ExpMap.count(pChr)==0)
        return _CHROM_();

    return (_CHROM_)ExpMap[pChr].getSubset(pStart,pEnd,options);
    }catch(...){throw;}
}

/** \brief Return a Chrom containing only the sites that overlap the given chr. Remove those elements from (this). REQUIRES COLLECTION TO BE SORTED
 *
* \param chr std::string : Name of scaffold to subset on
 * \param start int : Start position
 * \param end int : End position
 * \param options OverlapType : OverlapType if needed
 * \return _CHROM_ : Chrom containing the overlapping elements
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_CHROM_ uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::removeSubset(const std::string & pChr,const double pStart, const double pEnd, OverlapType overlap)
{
    try {
    if (ExpMap.count(pChr)==0)
        return _CHROM_();

    return (_CHROM_)ExpMap[pChr].removeSubset(pStart,pEnd,overlap);
    }catch(...){throw;}
}


/** \brief Return an EXP containing only the unarity sites that do not overlap does of the input structure. REQUIRES COLLECTION TO BE SORTED
 *
 *
 * \param compareExp uGenericNGSExperiment& Input,
 * \param options OverlapType How we determine if overlapping or not
 * \return uGenericNGSExperiment<_CHROM_,_BASE_> Experiment containing the sites in there appropriate Chroms
 *
 */

template<class _SELF_, typename _CHROM_, typename _BASE_>
_SELF_ uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::getDistinct(const std::string & pChr, const double pStart, const double pEnd,  OverlapType options)
{
    try{
    typename NGSExpMap::iterator iterMap;
    _SELF_ returnExp;
    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); iterMap++)
    {
        if (iterMap->first==pChr)
        {
            auto pChrom = this->getpChrom(iterMap->first);
            returnExp.addData(pChrom->getDistinct(pStart, pEnd) );
        }
        else
        {
            returnExp.addData(iterMap->second);
        }
    }
    return returnExp;
    }catch(...){throw;}
}

/** \brief Return an EXP containing only the unarity sites that do not overlap does of the input structure. This will also remove every chromosome that are not pChr. Remove the corresponding site form the Exp
 *
 *
 * \param compareExp uGenericNGSExperiment& Input,
 * \param options OverlapType How we determine if overlapping or not
 * \return uGenericNGSExperiment<_CHROM_,_BASE_> Experiment containing the sites in there appropriate Chroms
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_SELF_ uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::removeDistinct(const std::string & pChr, const double pStart, const double pEnd, OverlapType options)
{
    try{
    typename NGSExpMap::iterator iterMap;
    _SELF_ returnExp;
    _CHROM_ onlyChrom;
    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); iterMap++)
    {
        if (iterMap->first==pChr)
        {
            auto pChrom = this->getpChrom(iterMap->first);
            returnExp.addData(pChrom->removeDistinct(pStart, pEnd));
            onlyChrom=*pChrom;
        }
        else
        {
            returnExp.addData(iterMap->second);
        }
    }
    ExpMap.clear();
    if (onlyChrom.count())
        this->addData(onlyChrom);
    return returnExp;
    }catch(...){throw;}

}



/** \brief Receive a given chrom and add it to the Experiment
 *
 * \param inputChrom const _CHROM_& Chrom to add
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::addData(const _CHROM_ & inputChrom)
{
    /**< If chrom Already exist,  */
    try{
    if (ExpMap.count(inputChrom.getChr()) != 0)
    {
        _CHROM_* currentChrom;
        currentChrom=&(ExpMap[inputChrom.getChr()]);
        for (auto itChrom =inputChrom.begin(); itChrom!= inputChrom.end(); itChrom++)
        {
            currentChrom->addDataNoCheck(*itChrom);
        }
    }
    else /**< Make deep copy */
    {
        ExpMap.insert(std::pair<std::string,_CHROM_>(inputChrom.getChr(),inputChrom));
    }
    }catch(...){throw;}
}

/** \brief Receive a given chrom and add it to the Experiment. Replaces previous chr if any.
 *
 * \param inputChrom const _CHROM_& Chrom to add
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::replaceChr(const _CHROM_ & inputChrom)
{
    /**< Deep copy, create or replace as needed */
    ExpMap[inputChrom.getChr()]=inputChrom;
}

/** \brief Erase the structure associated with key, if the structure exists.
 *
 * \param pChrName const std::string& : Key of the scaffold to erase
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::removeChr(const std::string & pChrName)
{
    try {
	    uGenericNGSExperiment<_SELF_,_CHROM_, _BASE_>::ExpMap.erase(pChrName);
    }
    catch (...) {
        throw;
    }
}


/** \brief Return every element of THIS overlapping with parameter.
 *
 * \param uGenericNGSExperiment compareExp : Experiment to check for overlaps
 * \param type OverlapType What overlap to check
 * \return uGenericNGSExperiment<_CHROM_,_BASE_> Experiment containing all elements from THIS that overlap with parameter.
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_SELF_ uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::getOverlapping(_SELF_ &compareExp, OverlapType type)
{
    try {
    typename NGSExpMap::iterator iterMap;
    _CHROM_* pChrom;
    _SELF_ returnExp;
    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); iterMap++)
    {
        if (compareExp.isChrom(iterMap->first))
        {
            pChrom = compareExp.getpChrom(iterMap->first);
            returnExp.addData(iterMap->second.getOverlapping(*pChrom));
        }
    }
    return returnExp;
    }
    catch(...){throw;}
}

/** \brief Return every element of THIS overlapping with a specified chrom.
 *
 * \param uGenericNGSExperiment compareChrom : Experiment to check for overlaps
 * \param type OverlapType What overlap to check
 * \return uGenericNGSExperiment<_CHROM_,_BASE_> Experiment containing all elements from THIS that overlap with parameter.
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_SELF_ uGenericNGSExperiment<_SELF_, _CHROM_,_BASE_>::getOverlapping(_CHROM_ &compareChrom, OverlapType type)
{
    try
    {
        _SELF_ tempExp;
        tempExp.addData(compareChrom); // TODO: does not seem to work, but at least it compiles
        return getOverlapping(tempExp,type);
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Return every element of THIS overlapping with fixed interval.
 *
 * \param uGenericNGSExperiment compareChrom : Experiment to check for overlaps
 * \param type OverlapType What overlap to check
 * \return uGenericNGSExperiment<_CHROM_,_BASE_> Experiment containing all elements from THIS that overlap with parameter.
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
_SELF_ uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::getOverlapping(std::string chr, int start, int end, OverlapType type)
{
    try
    {
        _SELF_ tempChrom;
        auto item= _BASE_(chr,start,end);
        tempChrom.addData(item);
        return getOverlapping(tempChrom, type);
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Split each item into smaller equal size members and replace our vector of items with the new one.
 *
 * \param N int : Number of bins to make
 * \param type SplitType : How to manage splits that are not a multiple of N. Possible : STRICT - IGNORE - FILL - ADD
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::divideItemsIntoNBins(int N, SplitType type)
{
    try
    {
        for( auto& x : ExpMap)
        {
            x.second.divideItemsIntoNBins(N,type);
        }
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Failed in divideItemsIntoNBins"<<std::endl;
#endif
        throw uExp_operation_throw()<<string_error("Throwing while trying to call divideItemsIntoNBins() on all chroms");
    }
}

/** \brief Call divideItemsIntoBinofSize on every container item
 *
 * \param N int : Size of the members to build in BP
 * \param type SplitType : How to manage splits that are not a multiple of N. Possible : STRICT - IGNORE - FILL - ADD
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::divideItemsIntoBinofSize(int N, SplitType type)
{
    try
    {
        for(auto& x : ExpMap)
        {
            x.second.divideItemsIntoBinofSize(N,type);
        }
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Failed in divideItemsIntoBinofSize"<<std::endl;
#endif
        throw uExp_operation_throw()<<string_error("Throwing while trying to call divideItemsIntoBinofSize() on all chroms");
    }
}

/** \brief Remove sites from every chrom for which the predicate is true.
 *
 * \param pred UnaryPredicate Predicate to test, follows standard pattern
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryPredicate>
void uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::removeSpecificSites(UnaryPredicate pred)
{
    try
    {
        for(auto& x : ExpMap)
        {
            x.second.removeSpecificSites(pred);
        }
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Failed in divideItemsIntoBinofSize"<<std::endl;
#endif
        throw uExp_operation_throw()<<string_error("Throwing while trying to call removeSpecificSites() on all chroms");
    }

}

/** uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::\brief Accumulate information by querying all chromosomes
  *
  * Passes the given function to std::accumulate and runs it on every chrom element. Note that the passed
  * function needs to take directly a chrom structure as argument and an "accumulator" paramator. the function must return the new value
  * of the accumulator
  * \param binary_op BinaryOperation : Querying function to perform on the chromosomes collection
  * \param init InitialValue The initial value of the "accumulator"
  * \return The information accumulated by querying all the chromosomes
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class BinaryOperation, class InitialValue>
InitialValue uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::accumulateChromsInfo(BinaryOperation binary_op, InitialValue init) const
{
    // Force using sequential version for accumulate as parallel version
    // doesn't work if actual data type of InitialValue and _CHROM_ cannot be
    // converted back and forth.
   try{
    return __gnu_parallel::accumulate(std::begin(ExpMap), std::end(ExpMap), init,
                                      [&binary_op, &init] (decltype(init) partial, NGSExpPair element)
                                      -> decltype(binary_op(partial, element.second))
    {
        return binary_op(partial, element.second);
    }, __gnu_parallel::sequential_tag());
    }catch(...){throw;}
}

/** \brief Compute a value for all chromosomes in the experiment and return the resulting collection
  *
  * Takes the pased functions and passes it to std::transform. Stores the result in a Map structure
  * with each result mapping to the equivalent string. The passed function must return a non-void value.
  *
  * \param unary_op UnaryOperation : Unary operation to perform on all the chromosomes of the experiment
  * \return A collection of values computed on each chromosome by unary_op
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryOperation>
auto uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::computeOnAllChroms(UnaryOperation unary_op) const -> std::map<std::string, decltype(unary_op(_CHROM_()))>
{
    try {
    std::map<std::string, decltype(unary_op(_CHROM_()))> results;
    transform(std::begin(ExpMap), std::end(ExpMap), std::inserter(results, std::begin(results)), [&unary_op](NGSExpPair element)
    {
        return make_pair(element.first, unary_op(element.second));
    });
    return results;
    }catch(...){throw;}
}


/** \brief Transform the chromosomes collection by applying a certain function to all chromosomes
  *
  *  Takes the passed function and run them on every chromosome structure, via std::for_each.
  *
  * \param unary_op UnaryOperation : Unary operation to perform on the chromosomes collection
  * \return unary_op, the operation that was performed on all chromosomes
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryFunction>
UnaryFunction uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::applyOnAllChroms(UnaryFunction f)
{
    try{
    for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair& element)
    {
        f(element.second);
    });
    return f;
     }catch(...){throw;}
}

/** \brief Transform the chromosomes collection by applying a certain function to all chromosomes
  *
  * \param unary_op UnaryOperation : Unary operation to perform on the chromosomes collection
  * \return unary_op, the operation that was performed on all chromosomes
  * \sa applyOnAllChroms
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryFunction>
UnaryFunction uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::applyOnAllChroms(const UnaryFunction f)const
{
    try{
    for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair element)
    {
        f(element.second);
    });
    return f;
     }catch(...){throw;}
}

/** \brief Transform the sites of the EXP by applying a certain function
  *
  *  The passed function is given to applyOnAllSites for each chromosome structure.
  *
  * \param unary_op UnaryOperation : Unary operation to perform on the sites collection
  * \return unary_op, the operation that was performed on all sites
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryFunction>
UnaryFunction uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::applyOnSites(UnaryFunction f)
{
    try{
        for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair& element)
        {
            element.second.applyOnAllSites(f);
        });
        return f;
     }catch(...){throw;}
}

/**< Const version of its equivalent */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryFunction>
UnaryFunction uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::applyOnSites(const UnaryFunction f)const
{
    try{
    for_each(std::begin(ExpMap), std::end(ExpMap), [&f](const NGSExpPair& element)
    {
        element.second.applyOnAllSites(f);
    });
    return f;
     }catch(...){throw;}
}

/** \brief load data from Parser, convert to unitary and execute the given function. Does -not- necessarily add to EXP
 *
 * \param stream std::ifstream& file to load from
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryFunction>
void uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::loadWithParserAndRun(std::ifstream& pStream, std::string pType, UnaryFunction funct , int pBlockSize)
{
    try
    {
        uParser Curparser(&pStream, pType);
        loadWithParserAndRun(Curparser,funct,pBlockSize);
    }
    catch (...)
    {
        throw;
    }
}

template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryFunction>
/** \brief
 *
 * \param pParser uParser& Parser we should load from
 * \param funct UnaryFunction Function to run on each element
 * \param pBlockSize int Number of items to load at a time
 * \return void
 *
 */
void uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::loadWithParserAndRun(uParser& pParser, UnaryFunction funct , int pBlockSize)
{
    try
    {

        std::vector<uToken> loadedTokens;
        loadedTokens.resize(pBlockSize);
        while(!pParser.eof())
        {
            int curLoaded=0;
            /**< Load a block of data */
            while ((curLoaded<pBlockSize)&&(!pParser.eof()))
            {
                loadedTokens.at(curLoaded)=pParser.getNextEntry();

                curLoaded++;
            }
            /**< Operate */
            for(const uToken & curToken:loadedTokens)
            {
                funct( (_BASE_)(curToken) );
            }
        }
    }
    catch (...)
    {
        throw;
    }
}


/** \brief load data from Parser, convert to unitary and execute the given function. Does -not- necessarily add to EXP
 *
 * \param stream std::string filepath path to the file to load from
 * \return void
 *
 */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class UnaryFunction>
void uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::loadWithParserAndRun(std::string filepath, std::string pType, UnaryFunction f, int pBlockSize)
{
    try
    {
        uParser Curparser(filepath, pType);
        loadWithParserAndRun(Curparser,f,pBlockSize);
    }
    catch (...)
    {
        throw;
    }
}

/** \brief Count the chromosomes for which a certain predicate is true
  *
  * This function take a pointer to a predicate function; this function
  * pointer can either be a) * the name of a function taking a chromosome by
  * reference, b) a lambda function taking a chromosome by reference or c) a
  * member method of a chromosome using "mem_fun_ref". In all cases, the
  * function must return a boolean; true is the predicate is true, false
  * otherwise.
  *
  * \param p UnaryPredicate : Unary predicate to evaluate on all chromosomes
  * \return The number of chromosomes for which the predicate is true
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template <class UnaryPredicate>
typename std::iterator_traits<typename std::map<std::string,_CHROM_>::iterator>::difference_type
uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::countChromsWithProperty(UnaryPredicate pred) const
{
    try{
    return count_if(std::begin(ExpMap), std::end(ExpMap), [&pred](const NGSExpPair& element)
    {
        return pred(element.second);
    });
     }catch(...){throw;}
}

/** \brief Find the maximal chromosome according to a certain comparison
  *
  * This function take a pointer to a function to find the maximal chromosome;
  * this function pointer can either be a) the name of a function taking two
  * chromosomes as parameters, b) a lambda function taking two chromosomes as parameters
  * or c) a member method of a chromosome taking another chromosome as parameter using
  * "mem_fun_ref". In all cases, the function must return a boolean: true if
  * the first element is "lower" than the second, false otherwise.
  *
  * \param comp Compare : Binary comparison operation to perform on the chromosomes collection
  * \return An iterator to the maximal chromosome
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class Compare>
typename std::map<std::string,_CHROM_>::const_iterator uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::maxChrom(Compare comp) const
{
    return max_element(std::begin(ExpMap), std::end(ExpMap),
                       [&comp](const NGSExpPair& element1, const NGSExpPair& element2) -> bool
    {
        return comp(element1.second, element2.second);
    });
}

/** \brief Find the minimal chromosome according to a certain comparison
  *
  * This function take a pointer to a function to find the minimal chromosome;
  * this function pointer can either be a) the name of a function taking two
  * chromosomes as parameters, b) a lambda function taking two chromosomes as parameters
  * or c) a member method of a chromosome taking another chromosome as parameter using
  * "mem_fun_ref". In all cases, the function must return a boolean: true if
  * the first element is "lower" than the second, false otherwise.
  *
  * \param comp Compare : Binary comparison operation to perform on the chromosomes collection
  * \return An iterator to the minimal chromosome
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class Compare>
typename std::map<std::string,_CHROM_>::const_iterator uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::minChrom(Compare comp) const
{
    try{
    return min_element(std::begin(ExpMap), std::end(ExpMap),
                       [&comp](const NGSExpPair& element1, const NGSExpPair& element2) -> bool
    {
        return comp(element1.second, element2.second);
    });
     }catch(...){throw;}
}

/** \brief Find the minimal and maximal chromosomes according to a certain comparison
  *
  * This function take a pointer to a function to find the minimal and
  * maximal chromosomes; this function pointer can either be a) the name of a
  * function taking two chromosomes as parameters, b) a lambda function taking two
  * chromosomes as parameters or c) a member method of a chromosome taking another chromosome
  * as parameter using "mem_fun_ref". In all cases, the function must return
  * a boolean: true if the first element is "lower" than the second, false
  * otherwise.
  *
  * \param comp Compare : Binary comparison operation to perform on the chromosomes collection
  * \return A pair of iterators: the first indicates the minimal chromosome and the second, the maximal chromosome
  */
template<class _SELF_, typename _CHROM_, typename _BASE_>
template<class Compare>
std::pair<typename std::map<std::string,_CHROM_>::const_iterator, typename std::map<std::string,_CHROM_>::const_iterator> uGenericNGSExperiment<_SELF_,_CHROM_,_BASE_>::minAndMaxChroms(Compare comp) const
{
    try{
    return minmax_element(std::begin(ExpMap), std::end(ExpMap),
                          [&comp](const NGSExpPair& element1, const NGSExpPair& element2) -> bool
    {
        return comp(element1.second, element2.second);
    });
     }catch(...){throw;}
}

} // End of namespace NGS
#endif // UFORMATEXPERIMENT_H_INCLUDED
