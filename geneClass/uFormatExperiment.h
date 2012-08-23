#ifndef UFORMATEXPERIMENT_H_INCLUDED
#define UFORMATEXPERIMENT_H_INCLUDED
#include <fstream>
//#include "uParser.h"
//_BASE_ is our Tags, _CHROM_ our Chrom structure.
template<typename _CHROM_, typename _BASE_>
class uGenericNGSExperiment
{
    static_assert(
        std::is_convertible<_CHROM_, uGenericNGSChrom<_BASE_>>::value,
        "Both types do not use the same underlying data structures"
    );

    //We trust our template here is derived from uGenericNGSChrom
    typedef std::map<std::string,_CHROM_>      NGSExpMap;
    typedef typename NGSExpMap::iterator       NGSExpIter;
    typedef typename NGSExpMap::const_iterator NGSExpConstIter;
    typedef typename NGSExpMap::value_type     NGSExpPair;


protected:
    //Used for gradual mode



    int op_mode;
    std::ifstream* ourStream;
    std::map<std::string,_CHROM_>  ExpMap;
    void removeSite(std::string chr,int position);
    void inferChrSize();


public:

    virtual ~uGenericNGSExperiment(){};

    enum { DEFAULT = 0, GRADUAL = 1 };

    _CHROM_ getSubset(std::string chr, int start, int end, OverlapType options=OverlapType::OVERLAP_PARTIAL);
    _CHROM_ getDistinct(std::string chr, int start, int end, OverlapType options=OverlapType::OVERLAP_PARTIAL);

    //Should we be publicy allowing the turn of a pointer to our internal structure? I would assume not..


    const _CHROM_* getpChrom(const std::string & chrom) const
    {
        auto refer=&(ExpMap.find(chrom)->second);
        return (refer);
    };

    auto begin()const->decltype(ExpMap.cbegin()){return ExpMap.cbegin();};
    auto end()const->decltype(ExpMap.cend()){return ExpMap.cend();};
    auto begin()->decltype(ExpMap.begin()){return ExpMap.begin();};
    auto end()->decltype(ExpMap.end()){return ExpMap.end();};



    void combine(const _CHROM_ &);
    void addSite(const _BASE_ & newSite);
    long long count() const;

  //  int avgExpSiteSize();
  //  int minExpSiteSize();
  //  int maxExpSiteSize();

    int countExpUnique();

    //TODO Implement merge
    //void cutDuplicates();


    //TODO MAKE WORK
    bool isEndfile()
    {
        return ourStream->eof();
    };

    void setFileStream( std::ifstream& stream)
    {
        ourStream = &stream;
        op_mode = GRADUAL;
    };

    bool isModeGradual()
    {
        return op_mode == GRADUAL;
    };

    _BASE_ getSite(std::string chr, int position);
    int findPrecedingSite(std::string chr, int position, int low, int high);
    int getRegionCount(std::string chr, int start, int end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
    void sortData();


    virtual void loadFromTabFile(std::ifstream& stream);

    void writeAsBedFile(std::ostream& out);

    //TODO, code this whole part apparently? THIS IS WHAT HAPPENS WHEN WE DON'T TEST!
    template<typename _CHROMPAR_, typename _BASEPAR_>
    uGenericNGSExperiment getOverlapping(uGenericNGSExperiment<_CHROMPAR_,_BASEPAR_> &compareExp, OverlapType type=OverlapType::OVERLAP_COMPLETE);
    template<typename _BASEPAR_>
    uGenericNGSExperiment getOverlapping(uGenericNGSChrom<_BASEPAR_> &compareExp, OverlapType type=OverlapType::OVERLAP_COMPLETE);
    uGenericNGSExperiment getOverlapping(std::string chr, int start, int end, OverlapType type=OverlapType::OVERLAP_COMPLETE);


    uGenericNGSExperiment getDistinct(uGenericNGSExperiment &compareExp);
    /**<  ok from here*/


    void setChrSize(std::string chr, int chrSize){
        ExpMap[chr].setChromSize(chrSize);
    };
    int getChrSize(std::string chr)
    {
        return (ExpMap[chr].getChromSize());
    };







    int getSubsetCount(std::string chr, int start, int end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
    int getSubsetCount(uGenericNGS subsetReg, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);

    void divideItemsIntoBinofSize(int N, SplitType type=SplitType::STRICT);
    void divideItemsIntoNBins(int N, SplitType type=SplitType::STRICT);
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
    NGSExpMap getSpecificChroms(UnaryPredicate pred) const
    {
        NGSExpMap copyColl;
        copy_if(std::begin(ExpMap), std::end(ExpMap), std::inserter(copyColl, std::begin(copyColl)), [&pred](const NGSExpPair& element)
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
      * \param unary_op UnaryOperation : Unary operation to perform on the chromosomes collection
      * \return unary_op, the operation that was performed on all chromosomes
      */
    template<class UnaryFunction>
    UnaryFunction applyOnAllChroms(UnaryFunction f)
    {
        for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair element)
        {
            f(element.second);
        });
        return f;
    }

    template<class UnaryFunction>
    UnaryFunction applyOnSites(UnaryFunction f)
    {
        for_each(std::begin(ExpMap), std::end(ExpMap), [&f](NGSExpPair& element)
        {
            element.second.applyOnAllSites(f);
            //f(element.second);
        });
        return f;
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

  uGenericNGSExperiment():op_mode(DEFAULT) {};

    //TODO FIX THIS' MAKE EVERYTHIGN WITH IT PROTECTED
    _CHROM_* getpChrom(const std::string & chrom)
    {
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
    ptempChrom->addSite(newSite);
    }
    catch(std::exception & e)
    {
          #ifdef DEBUG
          cerr << "Catching and re-throwing in uFormatExp::addSite()"
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

/** \brief load basic data from file,
 *
 * \param stream std::ifstream& file to load from
 * \return void
 *
 *//*
template<typename _CHROM_, typename _BASE_>
 void uGenericNGSExperiment<_CHROM_, _BASE_>::loadFromFile(std::ifstream& stream, file_type p_fType ) */
//{

   /* uParser ourParse(streamm,p_fType);
  //  std::string tempString;
   // while(!std::getline(stream, tempString).eof())
  //  {
        while(!ourParse.eof())
        {
           addSite(Parser.getNextEntry());

        }

       addSite( static_cast<_BASE_>(factory::makeNGSfromTabString(tempString)));


   // }
} */

/** \brief Write our data as a legal bed file, filling only the first three columns
 *
 * \param out std::ofstream& stream to write to
 * \return void
 *
 */
template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_, _BASE_>::writeAsBedFile(std::ostream& out)
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
int uGenericNGSExperiment<_CHROM_,_BASE_>::getSubsetCount(std::string chr, int start, int end, OverlapType overlap)
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

//Return the number of elements in our experiment either overlap or included in the specified region.
template<typename _CHROM_, typename _BASE_>
int uGenericNGSExperiment<_CHROM_,_BASE_>::getSubsetCount(uGenericNGS subsetReg, OverlapType overlap)
{

    int count=0;
    count = getSubsetCount(subsetReg.getChr(),
                           subsetReg.getStart(),
                           subsetReg.getEnd());
    return count;
}


template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_, _BASE_>::sortData()
{
    //TODO Why won't this work with applyOnAllChroms
 for( auto it = this->begin(); it!=this->end(); it++)
            it->second.sortSites();;
 //  std::function<void (_CHROM_&)> funct;//=
 // funct=    (void(_CHROM_::*)()) &_CHROM_::sortSites;
  // applyOnAllChroms(funct);


  //  applyOnAllChroms(std::mem_fun_ref(static_cast<void (_CHROM_::*)()>(&_CHROM_::sortSites)));

}






template<typename _CHROM_, typename _BASE_>
int uGenericNGSExperiment<_CHROM_,_BASE_>::findPrecedingSite(std::string chr, int position, int low, int high)
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* tempChrom;

    tempChrom=&(ExpMap[chr]);

    return tempChrom->findPrecedingSite(position, low,high);
}

template<typename _CHROM_, typename _BASE_>
_BASE_ uGenericNGSExperiment<_CHROM_,_BASE_>::getSite(std::string chr, int position)
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* tempChrom;

    tempChrom=&(ExpMap[chr]);

    return tempChrom->getSite( chr,position);
}


//TODO make removeSubset and make a version that takes and reveices an NGSExp
template<typename _CHROM_, typename _BASE_>
_CHROM_ uGenericNGSExperiment<_CHROM_,_BASE_>::getSubset(std::string chr, int start, int end, OverlapType options)
{
    _CHROM_ returnChrom;
    _CHROM_* tempChrom;

    tempChrom=&(ExpMap[chr]);
    //If you want to use this, you will need to declare a constructur in the parent class of _CHROM_ to manage a _CHROM_<_BASE_> elementa
    //Copy constructor!
    returnChrom= (_CHROM_)tempChrom->getSubset(start,
                 end,
                 options);

    return returnChrom;
}

template<typename _CHROM_, typename _BASE_>
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getDistinct(uGenericNGSExperiment &compareExp)
{
    typename NGSExpMap::iterator iterMap;
    _CHROM_* pChrom;
    uGenericNGSExperiment<_CHROM_,_BASE_> returnExp;
    for (iterMap = ExpMap.begin(); iterMap != ExpMap.end(); iterMap++)
    {
        pChrom = compareExp.getpChrom(iterMap->first);
        returnExp.combine(iterMap->second.returnDistinct((uGenericNGSExperiment)compareExp));
    }
}


template<typename _CHROM_, typename _BASE_>
_CHROM_ uGenericNGSExperiment<_CHROM_,_BASE_>::getDistinct(std::string chr, int start, int end, OverlapType options)
{
    _CHROM_ returnChrom;
    _CHROM_* tempChrom;

    tempChrom=&(ExpMap[chr]);
    //If you want to use this, you will need to declare a constructur in the parent class of _CHROM_ to manage a _CHROM_<_BASE_> elementa
    //Copy constructor!
    returnChrom= (_CHROM_)tempChrom->getDistinct(start, end);

    return returnChrom;
}

template<typename _CHROM_, typename _BASE_>
void uGenericNGSExperiment<_CHROM_,_BASE_>::combine(const _CHROM_ & inputChrom)
{
    _CHROM_* currentChrom;
    std::vector<_BASE_> vecData;

    currentChrom=&(ExpMap[inputChrom.getChr()]);

    for(auto itChrom =inputChrom.begin(); itChrom!= inputChrom.end(); itChrom++){
        currentChrom->addSite(*itChrom);
    }
    /*
    vecData=inputChrom.returnVecData();

    for (unsigned int k=0; k < vecData.size(); k++)
    {
        tempChrom->addSite(vecData.at(k));
    } */
}

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
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getOverlapping(uGenericNGSChrom<_BASEPAR_> &compareChrom, OverlapType type)
{
    uGenericNGSExperiment<_CHROM_,_BASE_> tempExp;

    tempExp.combine(compareChrom);
    return getOverlapping(tempExp,type);

}

template<typename _CHROM_, typename _BASE_>
uGenericNGSExperiment<_CHROM_,_BASE_> uGenericNGSExperiment<_CHROM_,_BASE_>::getOverlapping(std::string chr, int start, int end, OverlapType type)
{
    uGenericNGSChrom<_BASE_> tempChrom;

    auto item= uGenericNGS(chr,start,end);
    auto castItem= static_cast<_BASE_>(item);
    tempChrom.addSite(castItem);
    return getOverlapping(tempChrom, type);
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
    catch(...)
    {
        std::cerr << "Failed in divideItemsIntoBinofSize"<<std::endl;
        throw;
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
    catch(...)
    {
        std::cerr << "Failed in divideItemsIntoBinofSize"<<std::endl;
        throw;
    }
}





namespace uGeneric{


    #define uGenericExp uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>,uGenericNGS>


}



#endif // UFORMATEXPERIMENT_H_INCLUDED
