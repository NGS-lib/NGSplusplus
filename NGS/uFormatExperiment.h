#ifndef UFORMATEXPERIMENT_H_INCLUDED
#define UFORMATEXPERIMENT_H_INCLUDED
#include <fstream>
#include "IO/Parser/uParser.h"

//

namespace NGS {
//_BASE_ is our Tags, _CHROM_ our Chrom structure.
enum class ReadMode{ DEFAULT, GRADUAL };
template<typename _CHROM_, typename _BASE_>
class uGenericNGSExperiment
{
    /**< Can be used to make sure EXP, chrom both use appropriate base */
    static_assert(
        std::is_convertible<_CHROM_, uGenericNGSChrom<_BASE_>>::value,
        "Both types do not use the same underlying data structures"
    );

    /**< Iterator typedef */
    typedef std::map<std::string,_CHROM_>      NGSExpMap;
    typedef typename NGSExpMap::iterator       NGSExpIter;
    typedef typename NGSExpMap::const_iterator NGSExpConstIter;
    typedef typename NGSExpMap::value_type     NGSExpPair;

    //TODO, const iterators public
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

    /**< Are we loading gradually? */
    ReadMode op_mode;
//    uParser m_parser;

    std::ifstream* ourStream;
    std::map<std::string,_CHROM_>  ExpMap;

    std::function<float(const _BASE_*)> sortGetStart=nullptr;
    std::function<float(const _BASE_*)> sortGetEnd=nullptr ;
    std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> m_comptFunc=compareStart;

    void removeSite(std::string chr,int position);

    void inferChrSize();


public:

    virtual ~uGenericNGSExperiment(){};
    uGenericNGSExperiment& operator=(const uGenericNGSExperiment& copFrom);
    uGenericNGSExperiment(const uGenericNGSExperiment&);

    auto begin()->decltype(ExpMap.begin()){return ExpMap.begin();};
    auto end()->decltype(ExpMap.end()){return ExpMap.end();};

    auto begin()const->decltype(ExpMap.cbegin()){return ExpMap.cbegin();};
    auto end()const->decltype(ExpMap.cend()){return ExpMap.cend();};

    _CHROM_ getSubset(std::string chr, float start, float end, OverlapType options=OverlapType::OVERLAP_PARTIAL);
   // _CHROM_ getDistinct(std::string chr, int start, int end, OverlapType options=OverlapType::OVERLAP_PARTIAL);

    //Should we be publicy allowing the turn of a pointer to our internal structure? I would assume not..

    _CHROM_ getChrom(const std::string & chrom) const
    {
        if (ExpMap.count(chrom)==0){
           throw ugene_operation_throw()<<string_error("Requested non-existent Chrom from Exp in getChrom(), value : " +chrom);
        }
        return  ExpMap.find(chrom)->second;;
    };


    const _CHROM_* getpChrom(const std::string & chrom) const
    {
        if (ExpMap.count(chrom)==0){
           throw ugene_operation_throw()<<string_error("Required pointer to non-existent Chrom from Exp in getpChrom(), value : " +chrom);
        }
        const auto refer=&(ExpMap.find(chrom)->second);
        return (refer);
    };

    void combine(const _CHROM_ &);
    void addSite(const _BASE_ & newSite);

    long long count() const;

    //Replace with parser
    bool isEndfile()
    {
        return ourStream->eof();
    };
    /**< For graduel loading */
    void setFileStream( std::ifstream& stream)
    {
        ourStream = &stream;
        op_mode = ReadMode::GRADUAL;
    };

    bool isModeGradual()
    {
        return op_mode == ReadMode::GRADUAL;
    };

    _BASE_ getSite(std::string chr, int position)const;
    _BASE_ getSite(typename std::vector<_BASE_>::const_iterator posItr)const;


    void sortSites();
    template<typename Compare>
    void sortSites(Compare comp,std::function<float(const _BASE_*)> getStart_funct=nullptr,std::function<float(const _BASE_*)> getEnd_funct=nullptr);
    bool isSorted()const;
    typename std::vector<_BASE_>::const_iterator findPrecedingSite(std::string chr, int position)const;
    typename std::vector<_BASE_>::const_iterator findNextSite(std::string chr, int position)const;

    virtual void loadFromTabFile(std::ifstream& stream);
    virtual void loadWithParser(std::ifstream&, std::string);
    virtual void loadWithParser(std::string, std::string);


    void writeAsBedFile(std::ostream& out)const;

    template<typename _CHROMPAR_, typename _BASEPAR_>
    uGenericNGSExperiment getOverlapping(uGenericNGSExperiment<_CHROMPAR_,_BASEPAR_> &compareExp, OverlapType type=OverlapType::OVERLAP_PARTIAL);
    template<typename _BASEPAR_>
    uGenericNGSExperiment getOverlapping(uGenericNGSChrom<_BASEPAR_> &compareExp, OverlapType type=OverlapType::OVERLAP_PARTIAL);
    uGenericNGSExperiment getOverlapping(std::string chr, int start, int end, OverlapType type=OverlapType::OVERLAP_PARTIAL);


    /**< Temp */
    void printChromSortStatus()const;

    uGenericNGSExperiment getDistinct(uGenericNGSExperiment &compareExp, OverlapType type=OverlapType::OVERLAP_PARTIAL);
    uGenericNGSExperiment getDistinct(_CHROM_ &compareChr, OverlapType type=OverlapType::OVERLAP_PARTIAL);
    uGenericNGSExperiment getDistinct(std::string chr, int start, int end, OverlapType type=OverlapType::OVERLAP_PARTIAL);
    uGenericNGSExperiment getDistinct(_BASE_ elem,  OverlapType options);
    /**<  ok from here*/


    void setChrSize(std::string chr, int chrSize){
        ExpMap[chr].setChromSize(chrSize);
    };
    int getChrSize(std::string chr)
    {
        return (ExpMap[chr].getChromSize());
    };

    int getSubsetCount(const std::string & chr, const float start, const float end, const OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
    int getSubsetCount(const uGenericNGS & subsetReg, const OverlapType overlap=OverlapType::OVERLAP_PARTIAL);

    void divideItemsIntoBinofSize(int N, SplitType type=SplitType::STRICT);
    void divideItemsIntoNBins(int N, SplitType type=SplitType::STRICT);


    /** \brief Remove sites from every chrom for which the predicate is true.
     *
     * \param pred UnaryPredicate Predicate to test, follows standard pattern
     * \return void
     *
     */
    template<class UnaryPredicate>
    void removeSpecificSites(UnaryPredicate pred){
    try
        {
        for(auto& x : ExpMap)
            {
                x.second.removeSpecificSites(pred);
            }
        }
        catch(uChrom_operation_throw &e)
        {
            #ifdef DEBUG
            std::cerr << "Failed in divideItemsIntoBinofSize"<<std::endl;
            #endif
            throw uExp_operation_throw()<<string_error("Throwing while trying to call removeSpecificSites() on all chroms");
        }

    }
    /**< Wrappers around the STL algorithms */
    /** \brief Accumulate information by querying all chromosomes
      *
      * This function take a pointer to a function to accumulate some information;
      * this function pointer can either be a) the name of a function taking two
      * parameters, an accumulator and a site or b) a lambda function taking two
      * parameters, an accumulator and a site. In all cases, the function must
      * return the new value of the accumulator.
      *
      * \param binary_op BinaryOperation : Querying function to perform on the chromosomes collection
      * \param init InitialValue The initial value of the "accumulator"
      * \return The information accumulated by querying all the chromosomes
      */
    template<class BinaryOperation, class InitialValue>
    InitialValue accumulateChromsInfo(BinaryOperation binary_op, InitialValue init) const
    {
        // Force using sequential version for accumulate as parallel version
        // doesn't work if actual data type of InitialValue and _CHROM_ cannot be
        // converted back and forth.
        return __gnu_parallel::accumulate(std::begin(ExpMap), std::end(ExpMap), init,
                                          [&binary_op, &init] (decltype(init) partial, NGSExpPair element)
                                          -> decltype(binary_op(partial, element.second))
        {
            return binary_op(partial, element.second);
        }, __gnu_parallel::sequential_tag());
    }

    /** \brief Compute a value for all chromosomes in the experiment and return the resulting collection
      *
      * This function take a pointer to a function to perform on all the
      * chromosomes in the collection; this function pointer can either be a)
      * the name of a function taking a chromosome by reference, b) a lambda
      * function taking a chromosome by reference or c) a member method of a
      * chromosome using "mem_fun_ref". In all cases, the function must return a
      * non void value.
      *
      * \param unary_op UnaryOperation : Unary operation to perform on all the chromosomes of the experiment
      * \return A collection of values computed on each chromosome by unary_op
      */
    template<class UnaryOperation>
    auto computeOnAllChroms(UnaryOperation unary_op) const -> std::map<std::string, decltype(unary_op(_CHROM_()))>
    {
        std::map<std::string, decltype(unary_op(_CHROM_()))> results;
        transform(std::begin(ExpMap), std::end(ExpMap), std::inserter(results, begin(results)), [&unary_op](NGSExpPair element)
        {
            return make_pair(element.first, unary_op(element.second));
        });
        return results;
    }

    /** \brief Compute a value for all chromosomes in the experiment and return the resulting collection
      *
      * This function take a pointer to a function to perform on all the
      * chromosomes in the collection; this function pointer can either be a)
      * the name of a function taking a chromosome by reference, b) a lambda
      * function taking a chromosome by reference or c) a member method of a
      * chromosome using "mem_fun_ref". In all cases, the function must return a
      * non void value.
      *
      * \param unary_op UnaryOperation : Unary operation to perform on all the chromosomes of the experiment
      * \return A collection of values computed on each chromosome by unary_op
      */
    template<class UnaryOperation>
    auto computeOnOneChrom(UnaryOperation unary_op, const std::string & chr) const -> std::map<std::string, decltype(unary_op(_CHROM_()))>
    {
        std::map<std::string, decltype(unary_op(_CHROM_()))> results;
       if (ExpMap.count(chr)){
            transform(std::begin(ExpMap), std::end(ExpMap), std::inserter(results, begin(results)), [&unary_op](NGSExpPair element)
            {
                return make_pair(element.first, unary_op(element.second));
            });
        }
        return results;
    }

    /** \brief Get the chromosomes for which a certain predicate is true
      *
      * This function take a pointer to a predicate function; this function
      * pointer can either be a) * the name of a function taking a chromosome by
      * reference, b) a lambda function taking a chromosome by reference or c) a
      * member method of a chromosome using "mem_fun_ref". In all cases, the
      * function must return a boolean; true is the predicate is true, false
      * otherwise.
      *
      * \param p UnaryPredicate : Unary predicate to evaluate on all chromosomes
      * \return A collection containing all the chromosomes for which the predicate is true
      */
    template<class UnaryPredicate>
    auto getSpecificChroms(UnaryPredicate pred) const->decltype(ExpMap)
   // NGSExpMap getSpecificChroms(UnaryPredicate pred) const
    {
//         auto begin()->decltype(ExpMap.begin()){return ExpMap.begin();};

       // NGSExpMap copyColl;
        decltype(ExpMap) copyColl;
        copy_if(std::begin(ExpMap), std::end(ExpMap), std::inserter(copyColl, std::begin(copyColl)), [&pred]( const typename decltype(ExpMap)::value_type& element)
        {
            return pred(element.second);
        });
        return copyColl;
    }

    /** \brief Transform the chromosomes collection by applying a certain function to all chromosomes
      *
      * This function take a pointer to a function to transform the chromosomes
      * collection; this function pointer can either be a) the name of a function
      * taking a site by reference, b) a lambda function taking a site by
      * reference or c) a member method of a site using "mem_fun_ref". In all
      * cases, the function must return void (any other return value will be
      * ignored).
      *
      *
      * \param unary_op UnaryOperation : Unary operation to perform on the chromosomes collection
      * \return unary_op, the operation that was performed on all chromosomes
      */
    template<class UnaryFunction>
    UnaryFunction applyOnAllChroms(UnaryFunction f)
    {
        for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair& element)
        {
            f(element.second);
        });
        return f;
    }

    /** \brief Transform the chromosomes collection by applying a certain function to one chromosomes
      *
      * This function take a pointer to a function to transform a single chrom
      ; this function pointer can either be a) the name of a function
      * taking a chrom by reference, b) a lambda function taking a chrom by
      * reference or c) a member method of a site using "mem_fun_ref". In all
      * cases, the function must return void (any other return value will be
      * ignored).
      *
      *
      * \param unary_op UnaryOperation : Unary operation to perform on the chromosomes collection
      * \return unary_op, the operation that was performed on all chromosomes
      */
    template<class UnaryFunction>
    UnaryFunction applyOnOneChrom(UnaryFunction f, const std::string & chr)
    {
        if (ExpMap.count(chr)){
            auto & elem = ExpMap[chr];
            f(elem);}
        return f;
    }

    /** \brief Transform the chromosomes collection by applying a certain function to all chromosomes
      *
      * This function take a pointer to a function to transform the chromosomes
      * collection; this function pointer can either be a) the name of a function
      * taking a chrom by reference, b) a lambda function taking a chrom by
      * reference or c) a member method of a chrom using "mem_fun_ref". In all
      * cases, the function must return void (any other return value will be
      * ignored).
      *
      *
      * \param unary_op UnaryOperation : Unary operation to perform on the chromosomes collection
      * \return unary_op, the operation that was performed on all chromosomes
      */
    template<class UnaryFunction>
    UnaryFunction applyOnAllChroms(const UnaryFunction f)const
    {
        for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair element)
        {
            f(element.second);
        });
        return f;
    }


    /** \brief Transform the sites of the EXP by applying a certain function
      *
      * This function take a pointer to a function to transform;
      * this function pointer can either be a) the name of a function
      * taking a site by reference, b) a lambda function taking a site by
      * reference or c) a member method of a site using "mem_fun_ref". In all
      * cases, the function must return void (any other return value will be
      * ignored).
      *
      * \param unary_op UnaryOperation : Unary operation to perform on the sites collection
      * \return unary_op, the operation that was performed on all sites
      */
    template<class UnaryFunction>
    UnaryFunction applyOnSites(UnaryFunction f)
    {
        for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair& element)
        {
            element.second.applyOnAllSites(f);
        });
        return f;
    }
     /**< Const version of its equivalent */
    template<class UnaryFunction>
    UnaryFunction applyOnSites(const UnaryFunction f)const
    {
        for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair& element)
        {
            element.second.applyOnAllSites(f);
            //f(element.second);
        });
        return f;
    }


/** \brief load data from Parser, convert to unitary and execute the given function. Does -not- necessarily add to EXP
 *
 * \param stream std::ifstream& file to load from
 * \return void
 *
 */
template<class UnaryFunction>
void loadWithParserAndRun(std::ifstream& pStream, std::string pType, UnaryFunction f , int pBlockSize=0)
{
    try {
    std::istream& refStream = pStream;
    uParser Curparser(&refStream, pType);

    while(!Curparser.eof()){
        int curLoaded=0;
        std::vector<uToken> loadedTokens;
            /**< Load a block of data */
            while ((curLoaded<pBlockSize)&&(!Curparser.eof()))
            {
                loadedTokens.push_back(Curparser.getNextEntry());
            }
            /**< Operate */
            for(const uToken & curToken:loadedTokens)
            {
                f( (_BASE_)(curToken) );
            }
        }
    }
    catch(...)
    { throw; }
}

template<class UnaryFunction>
 void loadWithParserAndRun(std::string filepath, std::string pType, UnaryFunction f  ,int pBlockSize=0)
{
    try {
    uParser Curparser(filepath, pType);
    while(!Curparser.eof()){
        int curLoaded=0;
        std::vector<uToken> loadedTokens;
            /**< Load a block of data */
            while ((curLoaded<pBlockSize)&&(!Curparser.eof()))
            {
                loadedTokens.push_back(Curparser.getNextEntry());
            }
            /**< Operate */
            for(const uToken & curToken:loadedTokens)
            {
                f( (_BASE_)(curToken) );
            }
        }
    }
    catch(...)
    { throw; }
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
    template <class UnaryPredicate>
    typename std::iterator_traits<NGSExpIter>::difference_type
    countChromsWithProperty(UnaryPredicate pred) const
    {
        return count_if(std::begin(ExpMap), std::end(ExpMap), [&pred](const NGSExpPair& element)
        {
            return pred(element.second);
        });
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
    template<class Compare>
    NGSExpConstIter maxChrom(Compare comp) const
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
    template<class Compare>
    NGSExpConstIter minChrom(Compare comp) const
    {
        return min_element(std::begin(ExpMap), std::end(ExpMap),
                           [&comp](const NGSExpPair& element1, const NGSExpPair& element2) -> bool
        {
            return comp(element1.second, element2.second);
        });
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
    template<class Compare>
    std::pair<NGSExpConstIter, NGSExpConstIter> minAndMaxChroms(Compare comp) const
    {
        return minmax_element(std::begin(ExpMap), std::end(ExpMap),
                              [&comp](const NGSExpPair& element1, const NGSExpPair& element2) -> bool
        {
            return comp(element1.second, element2.second);
        });
    }
    /**< End STL wrappers */

  uGenericNGSExperiment():op_mode(ReadMode::DEFAULT) {};

    //TODO Move to protect or private, should not be public
    _CHROM_* getpChrom(const std::string & chrom)
    {
        if (ExpMap.count(chrom)==0){
           throw ugene_operation_throw()<<string_error("Required pointer to non-existent Chrom from Exp in getpChrom(), value : " +chrom);
        }
        return &(ExpMap[chrom]);
    };

};

//Start uGenericNGSExperiment
template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_, _BASE_>::addSite(const _BASE_ & newSite)
{

    try {
    _CHROM_* ptempChrom;

    ptempChrom=&(ExpMap[newSite.getChr()]);
    ptempChrom->setChr(newSite.getChr());
    ptempChrom->addSite(newSite);
    }
    catch(std::exception & e)
    {
          #ifdef DEBUG
          std::cerr << "Catching and re-throwing in uFormatExp::addSite()" <<std::endl;
          #endif
        throw e;
    }
}

template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_, _BASE_>::inferChrSize(){
    applyOnAllChroms(std::mem_fun_ref(static_cast<void (_CHROM_::*)()>(&_CHROM_::inferChrSize)));
   //applyOnAllChroms(std::bind(static_cast<void (_CHROM_::*)()>(&_CHROM_::inferChrSize)));
}


/** \brief Remove a specific number from the specific subtype. //TODO Should this be public?
 *
 * \param chr std::string : Key to map element (chrom)
 * \param position int : position to remove
 * \return void uGenericNGSExperiment<_CHROM_,
 *
 */
template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_, _BASE_>::removeSite(std::string chr, int position)
{
    _CHROM_* tempChrom;
    tempChrom=&(ExpMap[chr]);
    try
    {
        if (position>=tempChrom.count())
        {

            abort();
        }
        tempChrom.removeSite(position);

    }
    catch(...)
    {
        std::cerr << "Crashing in removeSite() in uGenericNGSEXperiment,";
        throw;
    }
}

/** \brief load basic data from a tab delimited file, throw away the result.
 *
 * \param stream std::ifstream& file to load from
 * \return void
 *
 */
template<typename _CHROM_, typename _BASE_>
 void uGenericNGSExperiment<_CHROM_, _BASE_>::loadFromTabFile(std::ifstream& stream)
{
    std::string tempString;
    while(!std::getline(stream, tempString).eof())
    {
       addSite( static_cast<_BASE_>(factory::makeNGSfromTabString(tempString)));
    }
}

/** \brief load basic data from a Parser and load necessary data by passing to object constructor
 *
 * \param stream std::ifstream& file to load from
 * \return void
 *
 */
template<typename _CHROM_, typename _BASE_>
 void uGenericNGSExperiment<_CHROM_, _BASE_>::loadWithParser(std::ifstream& pStream, std::string pType)
{
    try {
    std::istream& refStream = pStream;
    uParser Curparser(&refStream, pType);
    while(!Curparser.eof()){
        addSite(Curparser.getNextEntry());
        }
    }
    catch(...)
    { throw; }
    inferChrSize();
}

template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_, _BASE_>::loadWithParser(std::string filepath, std::string pType)
{
    uParser ourParser(filepath, pType);
try {
        while (ourParser.eof()==false){
            this->addSite((ourParser.getNextEntry()));
    }
}
catch (...)
    {
        throw;
    }
    inferChrSize();
}




/** \brief Write our data as a legal bed file, filling only the first three columns
 *
 * \param out std::ofstream& stream to write to
 * \return void
 *
 */
template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_, _BASE_>::writeAsBedFile(std::ostream& out) const
{
    applyOnAllChroms(bind2nd(mem_fun_ref(&_CHROM_::outputBedFormat), out));
}

//Return the number of elements in our experiment
/** \brief Return our element count
 *
 * \return int : Number of elements
 *
 */
template<typename _CHROM_, typename _BASE_>
long long uGenericNGSExperiment<_CHROM_, _BASE_>::count() const
{
    return accumulateChromsInfo([](long long partialSum, _CHROM_ chrom) {return partialSum + chrom.count();}, 0LL);
}

/** \brief Return overlap of a specific position for a specific map
 *
 * \param chr std::string : Map Key ( chrom)
 * \param start int : Start position, must be positive
 * \param end int : End position, must be >= start
 * \param overlap int : Type of overlap
 * \return int
 *
 */
template<typename _CHROM_, typename _BASE_>
int uGenericNGSExperiment<_CHROM_,_BASE_>::getSubsetCount(const std::string & chr, const float start, const float end, OverlapType overlap)
{
    int count=0;
    typename NGSExpMap::iterator iterMap;
    _CHROM_* tempChrom;

    tempChrom=&(ExpMap[chr]);
    count = tempChrom->getSubsetCount(start,
                                      end,
                                      overlap);
    return count;
}

//TODO Deprecate or make general
template<typename _CHROM_, typename _BASE_>
int uGenericNGSExperiment<_CHROM_,_BASE_>::getSubsetCount(const uGenericNGS & subsetReg, const OverlapType overlap)
{

    int count=0;
    count = getSubsetCount(subsetReg.getChr(),
                           subsetReg.getStart(),
                           subsetReg.getEnd());
    return count;
}


/** \brief Sort Every site of every chrom based on location
 *
 * \return void uGenericNGSExperiment<_CHROM_,
 *
 */
template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_, _BASE_>::sortSites()
{
   sortSites(uGenericNGSExperiment<_CHROM_, _BASE_>::compareStart,&_BASE_::getStart,&_BASE_::getEnd);
}

  /** \brief Sort the sites vector by applying a certain comparison
      *      See documention on Chrom version for indepth comments
      *
      * \param comp Compare : Binary comparison operation to perform on the sites collection
      * \return void
      */
template<typename _CHROM_, typename _BASE_>
template<typename Compare>
void uGenericNGSExperiment<_CHROM_, _BASE_>::sortSites(Compare comp,std::function<float(const _BASE_*)> getStart_funct,std::function<float(const _BASE_*)> getEnd_funct)
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


template<typename _CHROM_, typename _BASE_>
/** \brief Returns false if at least one chrom is unsorted
 *
 * \return bool true if the experiment is sorted
 */
bool uGenericNGSExperiment<_CHROM_, _BASE_>::isSorted()const{
        bool sorted=true;
        applyOnAllChroms([&](_CHROM_& chrom){
                         if (chrom.isSorted(m_comptFunc)==false)
                         sorted=false;
                         });

    return sorted;
}

template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_,_BASE_>::printChromSortStatus()const
{
     for( auto it = this->begin(); it!=this->end(); it++){
            std::cerr << it->second.getChr() << " is sorted  " <<   it->second.getSortedStatus() <<std::endl;
    }
}

/** \brief Return an interator pointing to the element of the chr before or after the specified value
 *   Note that this is based on the current sort type so may not refer to genomic position.
 *   Requires the data to be sorted first
 * \param std::string chrom to search
 * \param int value to check
 * \return Iterator pointing to value or nullptr if invalid chr.
 *
 */
template<typename _CHROM_, typename _BASE_>
typename std::vector<_BASE_>::const_iterator uGenericNGSExperiment<_CHROM_,_BASE_>::findPrecedingSite(std::string chr, int position)const
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* tempChrom;
    if (ExpMap.count(chr))
        {
            tempChrom=&(ExpMap[chr]);
             return tempChrom->findPrecedingSite(position);
        }
return nullptr;
}
template<typename _CHROM_, typename _BASE_>
typename std::vector<_BASE_>::const_iterator uGenericNGSExperiment<_CHROM_,_BASE_>::findNextSite(std::string chr, int position)const
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* tempChrom;
     if (ExpMap.count(chr))
        {
            tempChrom=&(ExpMap[chr]);
            return tempChrom->findPrecedingSite(position);
        }
return nullptr;
}

/** \brief Get a specific site from a specific chrom. Overloaded to work with position or an interator, typically got from findPrecedingor findNext
 *
 * \param chr std::string
 * \param position int
 * \return _BASE_
 *
 */
template<typename _CHROM_, typename _BASE_>
_BASE_ uGenericNGSExperiment<_CHROM_,_BASE_>::getSite(std::string chr, int position)const
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* tempChrom;

    tempChrom=&(ExpMap[chr]);

    return tempChrom->getSite( chr,position);
}

template<typename _CHROM_, typename _BASE_>
_BASE_ uGenericNGSExperiment<_CHROM_,_BASE_>::getSite(typename std::vector<_BASE_>::const_iterator posItr)const
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* tempChrom;

    tempChrom=&(ExpMap[posItr->second.getChr()]);

    return tempChrom->getSite( posItr);
}






/** \brief Return a Chrom containing only the sites that overlap the given chr
 *
 *
 * \param chr std::string
 * \param start int
 * \param end int
 * \param options OverlapType
 * \return _CHROM_
 *
 */
 template<typename _CHROM_, typename _BASE_>
_CHROM_ uGenericNGSExperiment<_CHROM_,_BASE_>::getSubset(std::string chr, float start, float end, OverlapType options)
{
   if (ExpMap.count(chr)==0)
      return _CHROM_();

    return (_CHROM_)ExpMap[chr].getSubset(start,end,options);
}


/** \brief Return an EXP containing only the unarity sites that do not overlap does of the input structure
 *   Overloads include Experiment, Chrom, Site and chr,start,end
 *
 * \param compareExp uGenericNGSExperiment& Input, copies exist as mentionned prio
 * \param options OverlapType How we determine if overlapping or not
 * \return uGenericNGSExperiment<_CHROM_,_BASE_> Experiment containing the sites in there appropriate Chroms
 *
 */
 template<typename _CHROM_, typename _BASE_>
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getDistinct(uGenericNGSExperiment &compareExp, OverlapType options)
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* pChrom;
    uGenericNGSExperiment<_CHROM_,_BASE_> returnExp;
    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); iterMap++)
    {
        pChrom = compareExp.getpChrom(iterMap->first);
        returnExp.combine(iterMap->second.getDistinct((uGenericNGSExperiment)compareExp));
    }
    return returnExp;
}

template<typename _CHROM_, typename _BASE_>
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getDistinct(_CHROM_ &compareChr, OverlapType options)
{

    _CHROM_* pChrom;
    uGenericNGSExperiment<_CHROM_,_BASE_> returnExp;
    if (ExpMap.count(compareChr.getChr()) ){
        _CHROM_* tempChrom=&(ExpMap[compareChr.getChr()]);
        returnExp.combine(tempChrom->getDistinct(compareChr));
        }
    return returnExp;
}

template<typename _CHROM_, typename _BASE_>
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getDistinct(_BASE_ elem,  OverlapType options)
{
    typename NGSExpMap::iterator iterMap;
    uGenericNGSExperiment<_CHROM_,_BASE_> returnExp;
    _CHROM_* tempChrom=&(ExpMap[elem.getChr()]);

    return getDistinct(elem.getChr(),elem.getStart(),elem.getEnd(),options);
}

template<typename _CHROM_, typename _BASE_>
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getDistinct(std::string chr, int start, int end,  OverlapType options)
{
    typename NGSExpMap::iterator iterMap;
    uGenericNGSExperiment<_CHROM_,_BASE_> returnExp;
    _CHROM_* tempChrom=&(ExpMap[chr]);

    returnExp.combine(tempChrom->second.getDistinct(chr, start, end) );
    return returnExp;
}



/** \brief Receive a given chrom and add it to the Experiment
 *
 * \param inputChrom const _CHROM_& Chrom to add
 * \return void
 *
 */
template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_,_BASE_>::combine(const _CHROM_ & inputChrom)
{
    //TODO remove const ref to allow move semantics?
    _CHROM_* currentChrom;
    std::vector<_BASE_> vecData;

    currentChrom=&(ExpMap[inputChrom.getChr()]);
    for(auto itChrom =inputChrom.begin(); itChrom!= inputChrom.end(); itChrom++){
        currentChrom->addSite(*itChrom);
    }
}


/** \brief Return every element of THIS overlapping with parameter.
 *
 * \param uGenericNGSExperiment compareExp : Experiment to check for overlaps
 * \param type OverlapType What overlap to check
 * \return uGenericNGSExperiment<_CHROM_,_BASE_> Experiment containing all elements from THIS that overlap with parameter.
 *
 */
template<typename _CHROM_, typename _BASE_>
template<typename _CHROMPAR_, typename _BASEPAR_>
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getOverlapping(uGenericNGSExperiment<_CHROMPAR_,_BASEPAR_> &compareExp, OverlapType type)
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* pChrom;
    uGenericNGSExperiment<_CHROM_,_BASE_> returnExp;
    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); iterMap++)
    {
        pChrom = compareExp.getpChrom(iterMap->first);
        returnExp.combine(iterMap->second.getOverlapping(*pChrom));
    }
    return returnExp;
}

template<typename _CHROM_, typename _BASE_>
template<typename _BASEPAR_>
/** \brief Same as above, but compare a chrom rather then EXP
 *
 * \param uGenericNGSExperiment compareChrom : Experiment to check for overlaps
 * \param type OverlapType What overlap to check
 * \return uGenericNGSExperiment<_CHROM_,_BASE_> Experiment containing all elements from THIS that overlap with parameter.
 *
 */
 uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getOverlapping(uGenericNGSChrom<_BASEPAR_> &compareChrom, OverlapType type)
{
    try{
    uGenericNGSExperiment<_CHROM_,_BASE_> tempExp;

    tempExp.combine(compareChrom);
    return getOverlapping(tempExp,type);
    }
    catch(std::exception & e)
    {
        throw e;
    }
}
/**< See above, but using fixed positions */
template<typename _CHROM_, typename _BASE_>
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getOverlapping(std::string chr, int start, int end, OverlapType type)
{
    try {
    uGenericNGSChrom<_BASE_> tempChrom;

    auto item= uGenericNGS(chr,start,end);
    auto castItem= static_cast<_BASE_>(item);
    tempChrom.addSite(castItem);
    return getOverlapping(tempChrom, type);
    }
    catch(std::exception & e){
        throw e;
    }
}



/** \brief Split each item into smaller equal size members and replace our vector of items with the new one.
 *
 * \param N int : Number of bins to make
 * \param type SplitType : How to manage splits that are not a multiple of N. Possible : STRICT - IGNORE - FILL - ADD
 * \return void
 *
 */
template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_,_BASE_>::divideItemsIntoNBins(int N, SplitType type)
{
    try
    {
for( auto& x : ExpMap)
        {
            x.second.divideItemsIntoNBins(N,type);
        }
    }
    catch(uChrom_operation_throw &e)
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
template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_,_BASE_>::divideItemsIntoBinofSize(int N, SplitType type)
{
    try
    {
    for(auto& x : ExpMap)
        {
            x.second.divideItemsIntoBinofSize(N,type);
        }
    }
    catch(uChrom_operation_throw &e)
    {
        #ifdef DEBUG
        std::cerr << "Failed in divideItemsIntoBinofSize"<<std::endl;
        #endif
        throw uExp_operation_throw()<<string_error("Throwing while trying to call divideItemsIntoBinofSize() on all chroms");
    }
}



} // End of namespace NGS
#endif // UFORMATEXPERIMENT_H_INCLUDED
