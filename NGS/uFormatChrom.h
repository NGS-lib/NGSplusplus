#ifndef UFORMATCHROM_H_INCLUDED
#define UFORMATCHROM_H_INCLUDED
#include <climits>
#include <iostream>
#include <map>
#include <set>
#include "utility/utility.h"
#include <algorithm>
#include <parallel/numeric>
#include <functional>
#include <assert.h>
#include "uFormatBase.h"

/**********************************
      * The function pointer can either be a)
      * the name of a function taking a site by reference, b) a lambda
      * function taking a site by reference or c) a member method of a
      * function taking a site by reference or c) a member method of a
      * site using "mem_fun_ref" or "ref". In all cases, the function must return a
      * non void value.
      *
 ***************************/


namespace NGS
{
template<class _SELF_, typename _BASE_>

/******************
\brief The parent class for all Chrom structures
*
* This class is used as parent for all chromosome structures.
* It contains a vector of unatary base contig and the operations
* needed to manipulate them
*
********************/
class uGenericNGSChrom
{
    /**<  */
    static_assert(
        std::is_convertible<_BASE_, uGenericNGS<_BASE_>>::value,
        "The type does not inherit from uGenericNGS."
    );
    typedef std::vector<_BASE_> VecGenericNGS;
    typedef typename std::vector<_BASE_>::iterator VecGenIter;
    typedef typename std::vector<_BASE_>::const_iterator VecGenConstIter;


protected:

    std::vector<_BASE_> VecSites {}; /*!< std vector containing our unitary contigs  */
    std::string chr=""; /*!< Name of the scaffold/chromosome  */
    /**< Pointers to our functions and determines if sorted */
    bool m_isSorted=false; /*!< If we are in a sorted state or not */
    std::function<float(const _BASE_*)> sortGetStart=&_BASE_::getStart;  /*!<Pointer to the starting sort value */
    std::function<float(const _BASE_*)> sortGetEnd=&_BASE_::getEnd ;  /*!< Pointer to the end starting sort value. Typically set to equal Start */
    std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> m_comptFunc=compareStart; /*!< Pointer to sorting function */
    long long int chromSize=0; /*!< Size of the scaffold */

private :
    template <class Container>
    typename Container::iterator to_mutable_iterator(Container& c, typename Container::const_iterator it)
    {
        return c.begin() + (it - c.begin());
    }
protected:

    std::vector<long long> returnSiteSizes() const;
    /**< For internal use only, slightly faster */
    void addDataNoCheck( _BASE_ newSite);
    template<class L, typename S, typename R> friend class uGenericNGSExperiment;
    /**< Private site removal */
    void removeSite(int position);
    void removeSite(int start,int end);


public:

    _BASE_ getSite(long long pPos) const;
    virtual _SELF_ getCopy()const
    {
        assert (false);
    };
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


    /**< Constructors */
    uGenericNGSChrom() {};
    uGenericNGSChrom(const std::string & consString):chr(consString) {};
    uGenericNGSChrom(const std::string & consString, long int size);
    uGenericNGSChrom(const std::vector<_BASE_> &);
    uGenericNGSChrom(const std::string &, const std::vector<_BASE_> &);

    uGenericNGSChrom& operator=(const uGenericNGSChrom& copFrom);
    uGenericNGSChrom(const uGenericNGSChrom&);
    virtual ~uGenericNGSChrom<_SELF_,_BASE_> ()
    { ;  }

    /**< Utility wrappers to query collection */
    unsigned long long avgSiteSize() const;
    unsigned long long minSiteSize() const;
    unsigned long long maxSiteSize() const;
    unsigned long long sumSiteSize() const;
    void printStats(std::ostream& out) const;

    void inferChrSize();

    /**< In place, creates items */
    void divideItemsIntoNBins(int N, SplitType type=SplitType::STRICT);
    void divideItemsIntoBinofSize(int N, SplitType type=SplitType::STRICT);

    /**< Find according to sort value */
    VecGenConstIter findPrecedingSite(const float pPosition) const;
    VecGenConstIter findNextSite(const float pPosition) const;
    /**< removeSites overloads */
    void removeSite(VecGenConstIter pItrPos);
    void removeSite(VecGenConstIter pStartItr,VecGenConstIter pEndItr);
    /**< Functions to create and add items to our chrom */
    template <class _OTHER_>
    _BASE_ generateRandomSite(const int size, std::mt19937& engine, _OTHER_ exclList, const int sigma=0, const std::string ID="") const;
    _BASE_ generateRandomSite(const int size, std::mt19937& engine, const int sigma=0, const std::string ID="") const;
    template <class _OTHER_>
    void addNRandomSite(const int size, const int n, std::mt19937& engine, const _OTHER_ &exclList, const int sigma=0, const std::string ID="");
    void addNRandomSite(const int size, const int n, std::mt19937& engine, const int sigma=0, const std::string ID="");

    /**< Function to subset based on genomic positions and other collections */
    template <class _OTHER_>
    _SELF_ getOverlapping(_OTHER_ &pCompareChr,OverlapType pOverlap=OverlapType::OVERLAP_PARTIAL) const;
    template <class _OTHER_>
    long long int getOverlappingCount(_OTHER_ &pCompareChr,OverlapType pOverlap=OverlapType::OVERLAP_PARTIAL) const;

    template <class _OTHER_>
    _SELF_ getNotOverlapping(_OTHER_ &compareChr,OverlapType pOverlap=OverlapType::OVERLAP_PARTIAL) const;


    /**< Functions to manipulate generically ranges of our elements. Requires the collection to be sorted */
    _SELF_ getSubset(double p_start, double p_end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;
    long int getSubsetCount(double p_start, double p_end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL)const;

    _SELF_ removeSubset(double p_start, double p_end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
    _SELF_ removeDistinct(double p_start, double p_end, OverlapType options=OverlapType::OVERLAP_PARTIAL);

    _SELF_ getDistinct(double p_start, double p_end, OverlapType options=OverlapType::OVERLAP_PARTIAL) const;
    long int getDistinctCount(double p_start, double p_end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL)const;


    std::function<float(const _BASE_*)> getStartFunct() const;
    std::function<float(const _BASE_*)> getEndFunct() const;
    std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> getCompFunct() const;

    bool getSortedStatus() const;
    int count() const;

    long long int getChromSize() const;
    void setChromSize(long long int chromS);

    std::string getChr() const;
    void setChr(std::string pChr);

    std::vector<_BASE_> returnVecData() const;


    //TODO Vector overloads,
    // SELF_ getSubset(std::vector<std::pair<double,double>>, OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;
    // _SELF_ removeSubset(std::vector<std::pair<double,double>>, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
    //_SELF_ getDistinct(std::vector<std::pair<double,double>>, OverlapType options=OverlapType::OVERLAP_PARTIAL) const;
    //_SELF_ removeDistinct(std::vector<std::pair<double,double>>, OverlapType options=OverlapType::OVERLAP_PARTIAL);

    /**< Write using Writer */
    void writeWithWriter(uWriter& pWriter) const;

    /**< Function to add a unitary element */
    virtual void addData(const _BASE_ & newSite);
    virtual void addData(const uToken &);


    template<class UnaryOperation>
    std::vector<_BASE_> applyAndGetVecData(UnaryOperation unary_op);

    template<class UnaryOperation>
    _SELF_ applyAndGetChrom(UnaryOperation unary_op);


    template<class UnaryOperation>
    auto computeOnAllSites(UnaryOperation unary_op) -> std::vector<decltype(unary_op(_BASE_()))>;
    //TODO return chrom?,
    template<class UnaryPredicate>
    std::vector<_BASE_> getSpecificSites(UnaryPredicate pred) const;
    template<class UnaryPredicate>
    void removeSpecificSites(UnaryPredicate pred);

    template<class UnaryOperation>
    UnaryOperation applyOnAllSites(UnaryOperation f);

    template<class UnaryOperation>
    UnaryOperation applyOnAllSites(const UnaryOperation f) const;

    template<class BinaryOperation, class InitialValue>
    InitialValue accumulateSitesInfo(BinaryOperation binary_op, InitialValue init) const;
    /**< End STL wrappers */


    template<class Compare>

    void sortSites(Compare comp,std::function<float(const _BASE_*)> getStart_funct=nullptr,std::function<float(const _BASE_*)> getEnd_funct=nullptr);
    void sortSites();

    template<class Compare>
    bool isSorted(Compare comp) const;
    bool isSorted() const;


    long int countUnique() const;

    template<class Compare>
    VecGenConstIter minSite(Compare comp) const;
    template<class Compare>
    VecGenConstIter maxSite(Compare comp) const;
    template<class Compare>
    std::pair<VecGenConstIter, VecGenConstIter> minAndMaxSites(Compare comp) const;


    template <class UnaryPredicate>
    typename std::iterator_traits<VecGenIter>::difference_type countSitesWithProperty(UnaryPredicate p) const;
    /**< End STL wrappers */

    /**< Inline functions */
    /**< Public iterators */
    /** \brief Return an const iterator pointing to the first element of  VecSites
     *
     * \return Const random access iterator, pointing to the first element of VecSites
     *
     */
    auto begin()const->decltype(VecSites.cbegin())
    {
        return VecSites.cbegin();
    };
    /** \brief Return an const iterator pointing to the last element of  VecSites
     *
     * \return Const random access iterator, pointing to the last element of VecSites
     *
     */
    auto end()const->decltype(VecSites.cend())
    {
        return VecSites.cend();
    };
    /** \brief Return an iterator pointing to the first element of  VecSites
     *
     * \return Const random access iterator, pointing to the first element of VecSites
     *
     */
    auto begin()->decltype(VecSites.begin())
    {
        return VecSites.begin();
    };
    /** \brief Return an iterator pointing to the last element of  VecSites
     *
     * \return Const random access iterator, pointing to the last element of VecSites
     *
     */
    auto end()->decltype(VecSites.end())
    {
        return VecSites.end();
    };



};

/**< End inline functions */

/**<  Begin public */

/** \brief Construct with name and size
 *
 * \param consString std::string
 * \param size long int
 *
 */
template <class _SELF_, class _BASE_>
uGenericNGSChrom<_SELF_,_BASE_>::uGenericNGSChrom(const std::string & consString, long int size):chr(consString)
{
    try
    {
        setChromSize(size);
    }
    catch(std::exception & e)
    {
        throw e;
    }
}
template <class _SELF_, class _BASE_>
uGenericNGSChrom<_SELF_,_BASE_>::uGenericNGSChrom(const std::vector<_BASE_> & copyVec)
{
    for (_BASE_ elem: copyVec)
        addData(elem);
}

template <class _SELF_, class _BASE_>
uGenericNGSChrom<_SELF_,_BASE_>::uGenericNGSChrom(const std::string & pChrom, const std::vector<_BASE_> & pCopyVec):chr(pChrom)
{
    for (_BASE_ elem: pCopyVec)
        addData(elem);
}



template <class _SELF_,class _BASE_>
/** \brief Add a new element to our chrom
 *
 *  Add a new element to the collection. member value chr must  be the same for the collection and element.
 *
 * \exception ugene_exception_base : Thrown when the scaffold name of the element does not match the collection
 * \param newSite _BASE_ Ellement to add
 *
 */
void uGenericNGSChrom<_SELF_,_BASE_>::addData(const _BASE_& newSite)
{
    try
    {
        if (newSite.getChr()!=chr)
            throw ugene_exception_base()<<string_error("adding base to Chrom with non-matching scaffold/chr name");
        m_isSorted=false;
        VecSites.push_back(newSite);
    }
    catch(ugene_exception_base & e)
    {
        throw e;
    }
}
template <class _SELF_,class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::addData(const uToken & pToken)
{
    try
    {
        if (pToken.getParam(token_param::CHR)!=this->getChr())
            throw ugene_exception_base()<<string_error("adding base to Chrom with non-matching scaffold/chr name");
        m_isSorted=false;
        VecSites.push_back(_BASE_(pToken));
    }
    catch(ugene_exception_base & e)
    {
        throw e;
    }
}


/** \brief Write every element of the collection using the createToken() function and the passed Writer
 *
 * \param pWriter uWriter& The writer to use, will determine format.
 * \return void
 *
 */
template <class _SELF_,class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::writeWithWriter(uWriter& pWriter) const
{
    auto writeFunct= std::bind(&_BASE_::writeToOutput,std::placeholders::_1, pWriter);
    applyOnAllSites(writeFunct);
}

template <class _SELF_,class _BASE_>
/** \brief Overload for internal use, skips check
*
*
*   \exception std::exception : Thrown when
*/

void uGenericNGSChrom<_SELF_,_BASE_>::addDataNoCheck(_BASE_ newSite)
{
    m_isSorted=false;
    VecSites.push_back(std::move(newSite));
}

/** \brief Remove elements at a given position. Variants work similar, but erase a range or take iterators
 *
 * \param position int element to erase, starting from 0 to size
 * \return void
 *
 */
template <class _SELF_, class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::removeSite(int position)
{
    try
    {
        VecSites.erase(VecSites.begin()+(position));
    }
    catch(...)
    {
        throw;
    }

}
template <class _SELF_, class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::removeSite(int start, int end)
{
    try
    {
        m_isSorted=false;
        VecSites.erase(VecSites.begin()+(start),VecSites.begin()+(end));
    }
    catch(...)
    {
        throw;
    }
}




/** \brief Remove the range of elements pointed at by the iterator. This function is a relatively straightforward wrapping of
 *         std::erase and follows the same rules and canveats. Reminder, std:: erase will not erase the element pointed by the end-of-range iterator.
 *
 * \param start VecGenConstIter : Iterator pointing to the start of the range to erase
 * \param end VecGenConstIter : Iterator pointing to the end of the range to erase
 * \sa removeSite
 * \return void
 *
 */
template <class _SELF_, class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::removeSite(VecGenConstIter pStartItr,VecGenConstIter pEndItr)
{
    try
    {
        m_isSorted=false;
        /**< According to the C++11 standard, const iterator should be allowed in erase
        However, implementation does not seem to have caught up, hence this patch
         */

        VecSites.erase(to_mutable_iterator(VecSites,pStartItr),to_mutable_iterator(VecSites,pEndItr));
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Remove the element pointed at by iterator. This function is a relatively straightforward wrapping of
 *         std::erase and follows the same rules and canveats
 *
 * \param start VecGenConstIter : Iterator pointing to the start of the range to erase
 * \param end VecGenConstIter : Iterator pointing to the end of the range to erase
 * \sa removeSite
 * \return void
 *
 */
template <class _SELF_, class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::removeSite(VecGenConstIter pItrPos)
{
    try
    {
        m_isSorted=false;
        VecSites.erase(to_mutable_iterator(VecSites,pItrPos));
    }
    catch(...)
    {
        throw;
    }
}


/** \brief Generate a _BASE_ who's position (start) is randomized
 *
 *  Calling this function will generate a random start position between 0 and scaffoldSize
 *  and return a _BASE_ element using the position. Setting sigma larger then 0 will vary
 *  the size of the created element by sampling from a normal distribution.
 *
 *
 * \param size const int size of the region to generate
 * \param engine std::mt19937& type of engine to use, see manual for details
 * \param sigma const int Standard deviation if we wish to introduct randomness
 * \param ID const std::string Name of the string
 * \return _BASE_ the Element returned
 *
 */
template <class _SELF_,class _BASE_>
_BASE_ uGenericNGSChrom<_SELF_, _BASE_>::generateRandomSite(const int size, std::mt19937& engine, const int sigma, const std::string ID) const
{
    try
    {
        _SELF_ emptyExcl(this->getChr());
        return generateRandomSite(size,engine,emptyExcl,sigma,ID);
    }
    catch(param_throw &e)
    {
        throw e;
    }
}

template <class _SELF_,class _BASE_>
template <class _OTHER_>
/** \brief Generate a _BASE_ who's position (start) is randomized
 *
 * \param size const int size of the region to generate
 * \param engine std::mt19937& type of engine to use, see manual for details
 * \param excluList _OTHER_ : exclList is a container that implements getSubset. Items generated will not overlap any item in this container.
 * \param sigma const int Standard deviation if we wish to introduct randomness
 * \param ID const std::string Name of the string
 * \return _BASE_ the Element returned
 * \sa generateRandomSite
 */
_BASE_ uGenericNGSChrom<_SELF_,_BASE_>::generateRandomSite
(const int pSize,std::mt19937& engine,_OTHER_ exclList, const int sigma, const std::string ID) const
{
    try
    {
        _BASE_ returnTag;

        exclList.sortSites();

        bool found=false;
        int size = pSize;
        int max = this->getChromSize();
        if (size >=max)
            throw param_throw()<<string_error("Asked for element of size larger then scaffold in generateRandomSite()");



        while (!found)
        {
            {
                _BASE_ temptag;
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
                    temptag.setEnd(center+shift-1);
                    temptag.setStart(center-shift);
                    temptag.setChr(this->getChr());


                    if ((exclList.getSubset(temptag.getStart(),temptag.getEnd())).count()==0)
                    {
                        found=true;
                        returnTag=temptag;
                        returnTag.setChr(this->getChr());
                    }
                }
            }
        }

        return returnTag;
    }
    catch(param_throw & e)
    {
        throw e;
    }
}

template <class _SELF_,class _BASE_>
template <class T2>
/** \brief Will create N random elements and add them to the collection
 *
 * \param size const int :  size of the region to generate
 * \param N const int : Number of elements to generate
 * \param engine std::mt19937& : Type of engine to use, see manual for details
 * \param sigma const int : Standard deviation if we wish to introduct randomness
 * \param ID const std::string : Name of the string
 * \return _BASE_ : The Element returned
 * \sa generateRandomSite
 */
void uGenericNGSChrom<_SELF_,_BASE_>::addNRandomSite
(const int size, const int n, std::mt19937& engine, const T2& exclList, const int sigma, const std::string ID)
{
    try
    {
        std::string tempID=ID+"-"+getChr()+"-";
        for(int i=0; i<n; i++)
        {
            VecSites.push_back(generateRandomSite(size, engine,exclList,sigma, utility::concatStringInt(tempID, i)));
        }
        chr = VecSites.back().getChr();
    }
    catch(...)
    {
        throw;
    }

}

template <class _SELF_,class _BASE_>
/** \brief Will create N random elements and add them to the collection
*
* \param size const int :  size of the region to generate
* \param N const int : Number of elements to generate
* \param engine std::mt19937& : Type of engine to use, see manual for details
* \param sigma const int : Standard deviation if we wish to introduct randomness
* \param ID const std::string : Name of the string
* \return _BASE_ : The Element returned
* \sa addNRandomSite
*/
void uGenericNGSChrom<_SELF_,_BASE_>::addNRandomSite
(const int size, const int n, std::mt19937& engine, const int sigma, const std::string ID)
{
    try
    {
        _SELF_ emptyExcl;
        addNRandomSite(size,n, engine,emptyExcl,sigma, ID);
    }
    catch(...)
    {
        throw;
    }
}


/** \brief Return average element size
 *
 *  Wrapper function that will return the average size of the element sin the collection
 *
 *
 * \return unsigned long long : Average size of the elements
 *
 */
template <class _SELF_,class _BASE_>
unsigned long long uGenericNGSChrom<_SELF_,_BASE_>::avgSiteSize() const
{
    if (this->count() == 0)
        return 0;
    return sumSiteSize()/this->count();
}

template <class _SELF_, class _BASE_>
/** \brief return sum of sizes, including overlapping.
 *
 *  Wrapper function that calls accumulateSitesInfo and returns the sum
 *  of all elements lenghts.
 *
 * \return unsigned long long : Sum of the element lenghts.
 *
 */
unsigned long long uGenericNGSChrom<_SELF_,_BASE_>::sumSiteSize() const
{
    return accumulateSitesInfo([](unsigned long long partialSum, _BASE_ item) -> unsigned long long
    {
        return partialSum + item.getLenght();
    }, 0ULL);
}

template <class _SELF_, class _BASE_>
/** \brief Returns size of smallest element
 *
 *  Wrapper function that returns the size of the smallest element in
 *  the collection.
 *
 *
 * \return unsigned long long : Size of the smallest element.
 */
unsigned long long uGenericNGSChrom<_SELF_ ,_BASE_>::minSiteSize() const
{
    if (this->count() == 0)
        return 0;
    return minSite(compareLenght)->getLenght();
}

template <class _SELF_, class _BASE_>

/** \brief Return the size of the largest element
 *
 *  Wrapper function that returns the size of the largest element in the collection.
 *
 * \return unsigned long long: Size of the largest element.
 *
 */
unsigned long long uGenericNGSChrom<_SELF_,_BASE_>::maxSiteSize() const
{
    if (this->count() == 0)
        return 0;
    return maxSite(compareLenght)->getLenght();
}


template <class _SELF_, class _BASE_>
/** \brief Return how unique start/end combinations are present in the collection
 *
 *  The function will return the number of unique start/end combinations. Ex: A collection contains two
 *  elements at position 100-200 and three elements at other positions. The function will return 4 as the two
 *  elements at position 100-200 will be counted once.
 *
 * \return long int Number of unique start/end positions.
 *
 */
long int uGenericNGSChrom<_SELF_,_BASE_>::countUnique() const
{
    std::pair<long long int,long long int> current;
    std::set< std::pair<long long int,long long int> >  UniqueSet;

    for (auto iterVec = VecSites.begin() ; iterVec!= VecSites.end(); iterVec++)
    {
        current.first= iterVec->getStart();
        current.second = iterVec->getEnd();
        UniqueSet.insert(current);
    }
    return UniqueSet.size();

}

template <class _SELF_,class _BASE_>
/** \brief Returns a vector that contains the size of every element in the collection.
 *
 * \return std::vector<long long>: Vector containing the lenght of every element.
 *
 */
std::vector<long long> uGenericNGSChrom<_SELF_,_BASE_>::returnSiteSizes() const
{
    return computeOnAllSites([] (_BASE_ elem) -> long long {return elem.getLenght();});
}


template <class _SELF_,class _BASE_>
/** \brief Wrapper class, outputs various information on the collection to given stream
 *
 * \param out std::ostream& Stream to output to, can be standard output.
 * \return void
 *
 */
void uGenericNGSChrom<_SELF_,_BASE_>::printStats(std::ostream& out) const
{
    typename std::vector<long long> quarts;
    /**< Get a vector containing the lenght of every site */
    quarts = utility::quartilesofVector(computeOnAllSites([] (_BASE_ elem) -> long long {return elem.getLenght();}));

    out <<"Number of sites"<< "\t"<< this->count()<<"\n";
    out <<"Average sites size:"<< "\t"<< this->avgSiteSize()<<"\n";
    out <<"Median size: "<< "\t"<< quarts.at(1)<<"\n";
    out <<"q1 :" << "\t"<< quarts.at(0)<<"\n";
    out <<"q3 :" << "\t"<< quarts.at(2)<<"\n";
    auto minAndMax = minAndMaxSites(compareLenght);
    out <<"Min sites size:"<< "\t"<< minAndMax.first->getLenght() <<"\n";
    out <<"Max sites size:"<< "\t"<< minAndMax.second->getLenght() <<"\n";

}


/** \brief Find first item preceding a given value for our current sort type. Data must be sorted
 *
 *  This function will find the position preceding the given value as valid to the current sort condition. By default,
 *  the sort condition is based on the start position of every element, but proper use of the sorting function can modify this condition.
 *
 *  If called on an unsorted collection, this function will raise an exception.
 *
 * \exception unsorted_throw : Will throw if the collection is unsorted
 * \exception ugene_exception_base : Will throw if the proper getters where not set
 * \param position const int& value to evaluate from (based on the sort type: i.e.: if sorted by score, it will use score as position)
 * \return typename std::vector<_BASE_>::const_iterator Constant iterator pointing to the item
 *
 */
template <class _SELF_,class _BASE_>
typename std::vector<_BASE_>::const_iterator uGenericNGSChrom<_SELF_,_BASE_>::findPrecedingSite(const float pPosition) const
{
    try
    {
        /**< If unsorted, fail */

        if (VecSites.size()==0)
            return VecSites.end();
        if (m_isSorted==false)
            throw unsorted_throw() <<string_error("findPrecedingSite called on unsorted vector \n") ;
        if ((sortGetStart==nullptr)||(sortGetEnd==nullptr))
            throw ugene_exception_base() <<string_error(" findPrecedingSite called on chrom without appropriate start or end function\n") ;
        auto comp = [&] (const _BASE_ &item1, const float &item2)
        {
            return sortGetStart(&item1)< item2;
        };

        /**< Compare, sort Value */
        auto lower = std::lower_bound(VecSites.begin(), VecSites.end(), pPosition, comp);

        /**< If result is our first item, then no item precedes it */
        if ((lower==VecSites.begin()))
            return VecSites.end();
        /**< If result is end, every idem precedes the value */
        /**<Return item precedes and as such is LESS then pPosition. if no item was found, last item is closest to value  */
        lower--;
        return (lower);
    }
    catch (unsorted_throw & e)
    {
#ifdef DEBUG
        std::cerr << "FindPrecedingSite called on unsorted vector" <<std::endl;
#endif
        throw e;
    }
    catch (ugene_exception_base & e)
    {
#ifdef DEBUG
        std::cerr << "Calling findPrecedingSite and you did not provide an aproriate get function" <<std::endl;
#endif
        throw e;
    }
}

/** \brief Find first element after a given value for our current sort type. Data must be sorted
 *
 *  This function will find the element after the given value as valid to the current sort condition. By default,
 *  the sort condition is based on the start position of every element, but proper use of the sorting function can modify this condition.
 *
 *  If called on an unsorted collection, this function will raise an exception.
 *
 * \exception unsorted_throw : Will throw if the collection is unsorted
 * \exception ugene_exception_base : Will throw if the proper getters where not set
 * \param position const int& value to evaluate from (based on the sort type: i.e.: if sorted by score, it will use score as position)
 * \return typename std::vector<_BASE_>::const_iterator Constant iterator pointing to the value
 *
 */
template <class _SELF_,class _BASE_>
typename std::vector<_BASE_>::const_iterator uGenericNGSChrom<_SELF_,_BASE_>::findNextSite(const float pPosition) const
{
    try
    {
        /**< If unsorted, fail */

        if (VecSites.size()==0)
            return VecSites.end();

        if (m_isSorted==false)
            throw unsorted_throw() <<string_error("findNextSite called on unsorted vector \n") ;
        if ((sortGetStart==nullptr)||(sortGetEnd==nullptr))
            throw ugene_exception_base() <<string_error(" findNextSite called on chrom without appropriate start or end function\n") ;

        /**< Return true comparitor if item1 smaller then item 2 */
        auto comp = [&] (const float &item1, const _BASE_ &item2)
        {
            return item1< sortGetStart(&item2);
        };

        /**< Compare, sort Value */
        auto upper = std::upper_bound(VecSites.begin(), VecSites.end(), pPosition, comp);

        /**< If no result*/
        if (upper==VecSites.end())
            return VecSites.end();

        /**<Return the item greater then value*/
        return (upper);
    }
    catch (unsorted_throw & e)
    {
#ifdef DEBUG
        std::cerr << "findNextSite called on unsorted vector" <<std::endl;
#endif
        throw e;
    }
    catch (ugene_exception_base & e)
    {
#ifdef DEBUG
        std::cerr << "Calling findNextSite and you did not provide an aproriate get function" <<std::endl;
#endif
        throw e;
    }
}

/** \brief Return the count of a data subset, based on the current sort type.
 *
 *  Will return the number of elements in the specified interval. This function requires the elements to be sorted
 *  and by default the sort condition is based on genomic position ( start/end ).
 *
 *
 * \param start double start interval value
 * \param end double end interval value
 * \param overlap OverlapType OverlapType from Enum
 * \return int Number of elements in our range
 *
 */
template <class _SELF_,class _BASE_>
long int uGenericNGSChrom<_SELF_,_BASE_>::getSubsetCount(double p_start, double p_end, OverlapType overlap) const
{
    try
    {
        auto pos = this->findPrecedingSite(p_start);
        /**<  If no tag leftwise, we start at beginning*/
        if (pos==this->end())
            pos=this->begin();
        int tagcount=0;
        for (; pos != this->end(); pos++)
        {
            if (sortGetStart(&(*pos))> p_end)
                break;
            if (utility::isOverlap(sortGetStart(&(*pos)), sortGetEnd(&(*pos)),p_start, p_end,overlap))
                tagcount++;
        }

        return tagcount;
    }
    catch (unsorted_throw & e)
    {
        throw e;
    }
    catch (ugene_exception_base & e)
    {
        throw e;
    }
}


/** \brief Return a subset of our data that overlaps range start/end, based on current sort type.
 *
 * \param start double: Start of interval
 * \param end double: End of interval
 * \param overlap OverlapType: Type of overlap
 * \return uGenericNGSChrom<_BASE_>: Chrom structure containing our element subset
 *
 */
template <class _SELF_,class _BASE_>
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::getSubset(double p_start, double p_end, OverlapType overlap) const
{
    try
    {
        _SELF_ returnChrom;
        returnChrom.setChr(this->getChr());

        auto pos = this->findPrecedingSite(p_start);

        /**<  If no tag leftwise, we start at beginning*/
        if (pos==this->end())
            pos=this->begin();

        for (; pos != this->end(); pos++)
        {
            if (sortGetStart(&(*pos))> p_end)
                break;
            _BASE_ temp;
            if (utility::isOverlap(sortGetStart(&(*pos)), sortGetEnd(&(*pos)),p_start, p_end,overlap))
                returnChrom.addDataNoCheck(*pos);
        }
        return returnChrom;
    }
    catch (unsorted_throw & e)
    {
        throw e;
    }
    catch (ugene_exception_base & e)
    {
        throw e;
    }
}

/** \brief Remove a subset of elements based on the current comparison value. Returns the removed elements as a new container
 *
 * \param p_start float : Start of interval
 * \param p_end float : End of internval
 * \param overlap OverlapType : Type of comparison
 * \return uGenericNGSChrom<_BASE_> : Chrom structure containing the removed elements
 *
 */
template <class _SELF_,class _BASE_>
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::removeSubset(double p_start, double p_end, OverlapType overlap)
{
    try
    {
        _SELF_ returnChrom;
        returnChrom.setChr(this->getChr());

        auto posIter = this->findPrecedingSite(p_start);

        bool start=true;
        auto startPosIter=posIter;
        auto endPosIter= posIter;
        /**<  If no tag leftwise, we start at beginning*/
        if (posIter==this->end())
            posIter=this->begin();

        for (; posIter != this->end(); posIter++)
        {
            if (sortGetStart(&(*posIter))> p_end)
                break;

            if (utility::isOverlap(sortGetStart(&(*posIter)), sortGetEnd(&(*posIter)),p_start, p_end,overlap))
            {
                if (start)
                {
                    start =false;
                    startPosIter=posIter;
                    endPosIter= posIter;
                    returnChrom.addDataNoCheck(*posIter);
                }
                else
                {
                    endPosIter=posIter;
                    returnChrom.addDataNoCheck(*posIter);
                }
            }
        }
        if(!start)
        {

            endPosIter++;
            /**< std erase format, so endPos is not erased. */
            this->removeSite(startPosIter,endPosIter);
        }

        return returnChrom;
    }
    catch (unsorted_throw & e)
    {
        throw e;
    }
    catch (ugene_exception_base & e)
    {
        throw e;
    }
}

/** \brief Remove all elements not in an interval based on the current comparison value. Returns the removed elements as a new container
 *
 * \param p_start float : Start of interval
 * \param p_end float : End of internval
 * \param overlap OverlapType : Type of comparison
 * \return uGenericNGSChrom<_BASE_> : Chrom structure containing the removed elements
 *
 */
template <class _SELF_,class _BASE_>
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::removeDistinct(double p_start, double p_end, OverlapType options)
{
    try
    {

        _SELF_ returnChrom;
        returnChrom.setChr(this->getChr());

        auto posIter = this->findPrecedingSite(p_start);
        /**< Backtrack if needed */

        while (posIter!=this->begin() && sortGetEnd(&(*posIter))>p_start )
        {
            posIter--;
        }

        /**<  Remove elements before our interval, if any*/
        if (posIter!=this->end())
        {
            for (auto curCopy=this->begin(); curCopy!=posIter; curCopy++)
                returnChrom.addDataNoCheck(*curCopy);

            /**< Validate the element itself */
            if (!(utility::checkOverlap(sortGetStart(&(*posIter)),sortGetEnd(&(*posIter)), p_start, p_end)))
            {
                returnChrom.addDataNoCheck(*posIter);
                posIter++;
                std::cout<< "removed border"<<std::endl;
            }
            std::cout<< "removed before"<<std::endl;
            this->removeSite(this->begin(),posIter);
        }

        /**< Check elements until past the interval */
        for (posIter=this->begin(); posIter != this->end(); posIter++)
        {
            if (sortGetStart(&(*posIter))> p_end)
                break;
        }
        /**< Past the elements we keep, remove all */
        if(posIter!=this->end())
        {
            for (auto curCopy=posIter; curCopy!=this->end(); curCopy++)
                returnChrom.addDataNoCheck(*curCopy);
            std::cout<< "removed after"<<std::endl;
            this->removeSite(posIter,this->end());
        }
        return returnChrom;
    }
    catch (unsorted_throw & e)
    {
        throw e;
    }
    catch (ugene_exception_base & e)
    {
        throw e;
    }

}

/**< Return elements of A that overlap B */
/** \brief Wrapper function that returns a chrom structure containing the elements of that overlap another chrom structur
 *
 *      This function return a collection. This collection contains every element of THIS that overlaps an element of pCompareChr. This comparison
 *      is always based on genomic positions ( start end )
 *
 * \param _OTHER_ & pCompareChr : A compatible chrom collection
 * \param OverlapType overlap  : type of overlap
 * \return Chrom collection containing all overlapping elements.
 * \sa getNotOverlapping
 * \sa getOverlappingCount
 */
template <class _SELF_,class _BASE_>
template <class _OTHER_>
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::getOverlapping(_OTHER_ &pCompareChr,OverlapType pOverlap) const
{
    _SELF_ returnChr;

    if (getChr()!=pCompareChr.getChr())
        return returnChr;

    for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
    {
        for(auto compit= pCompareChr.begin(); compit!=pCompareChr.end(); compit++)
        {
            if (utility::isOverlap(it->getStart(), it->getEnd(),compit->getStart(),compit->getEnd(),pOverlap))
            {
                returnChr.addDataNoCheck(*it);
                break;
            }
        }
    }
    return returnChr;
}


/**< Return the number of elements of A that overlap B */
/** \brief Wrapper function that returns the number of elements of that overlap another chrom structur
 *
 *      This function return a collection. This collection contains every element of THIS that overlaps an element of pCompareChr. This comparison
 *      is always based on genomic positions ( start end )
 *
 * \param _OTHER_ & pCompareChr : A compatible chrom collection
 * \param OverlapType overlap  : type of overlap
 * \return Chrom collection containing all overlapping elements.
 * \sa getOverlapping
 */
template <class _SELF_,class _BASE_>
template <class _OTHER_>
long long int uGenericNGSChrom<_SELF_,_BASE_>::getOverlappingCount(_OTHER_ &pCompareChr,OverlapType pOverlap) const
{

    long long int overlapCount=0;
    if (getChr()!=pCompareChr.getChr())
        return overlapCount;

    for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
    {
        for(auto compit= pCompareChr.begin(); compit!=pCompareChr.end(); compit++)
        {
            if (utility::isOverlap(it->getStart(), it->getEnd(),compit->getStart(),compit->getEnd(),pOverlap))
            {
                overlapCount++;
                break;
            }
        }
    }
    return overlapCount;
}



/**< Return elements of A that overlap B */
/** \brief Wrapper function that returns a chrom structure containing the elements that do notoverlap another chrom structur
 *
 *      This function return a collection. This collection contains every element of THIS that overlaps an element of compareChr. This comparison
 *      is always based on genomic positions ( start end )
 *
 * \param _OTHER_ &compareChr : An compatible chrom collection
 * \param overlap : type of overlap
 * \return Chrom collection containing all overlapping elements.
 * \sa getOverlapping
 */
template <class _SELF_,class _BASE_>
template <class _OTHER_>
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::getNotOverlapping(_OTHER_ &pCompareChr,OverlapType pOverlap)const
{
    try
    {

        if (getChr()!=pCompareChr.getChr())
        {
            return this->getCopy();
//            return (returnChr=*this);
        }
        _SELF_ returnChr;
        bool add=true;
        for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
        {
            for(auto compit= pCompareChr.begin(); compit!=pCompareChr.end(); compit++)
            {
                if (utility::isOverlap(it->getStart(), it->getEnd(),compit->getStart(),compit->getEnd(),pOverlap))
                {
                    add=false;
                    break;
                }
            }
            if (add)
                returnChr.addDataNoCheck(*it);
            add =true;
        }
        return returnChr;
    }
    catch(std::exception & e)
    {
        throw e;
    }
}

/** \brief Removes the elements that overlap the given interval and return a Chrom containing the leftover. Based on the current sort.
 *
 * \param p_start double starting value of interval
 * \param p_end double end value of the interval
 * \param overlap OverlapType How to determine if the overlap is true or false
 * \return _SELF_ Chrom without the elements that overlap the interval
 */
template <class _SELF_,class _BASE_>
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::getDistinct(double p_start, double p_end, OverlapType overlap) const
{
    try
    {
        _SELF_ returnChrom;
        returnChrom= this->getCopy();
        returnChrom.sortGetStart=sortGetStart;
        returnChrom.sortGetEnd=sortGetEnd;
        auto posIter = returnChrom.findPrecedingSite(p_start);

        /**<  If no tag leftwise, we start at beginning*/
        if (posIter==returnChrom.end())
            posIter=returnChrom.begin();
        bool start=true;
        auto startPosIter=posIter;
        auto endPosIter= posIter;
        for (; posIter != returnChrom.end(); posIter++)
        {
            if (sortGetStart(&(*posIter))> p_end)
                break;

            if (utility::isOverlap(sortGetStart(&(*posIter)), sortGetEnd(&(*posIter)),p_start, p_end,overlap))
            {
                if (start)
                {
                    start =false;
                    startPosIter=posIter;
                    endPosIter= posIter;
                }
                else
                {
                    endPosIter=posIter;
                }
            }
        }
        if(!start)
        {
            /**< std erase format, so endPos is not erased. */
            endPosIter++;
            returnChrom.removeSite(startPosIter,endPosIter);
        }

        return returnChrom;
    }
    catch (unsorted_throw & e)
    {
#ifdef DEBUG
        std::cerr << "getDistinct called on unsorted vector" <<std::endl;
#endif
        throw e;
    }
    catch (ugene_exception_base & e)
    {
#ifdef DEBUG
        std::cerr << "Calling getDistinct and you did not provide an aproriate get function" <<std::endl;
#endif
        throw e;
    }
}


/** \brief Count the number of element that do not overlap the given interval based on the current sort.
 *
 * \param p_start double starting value of interval
 * \param p_end double end value of the interval
 * \param overlap OverlapType How to determine if the overlap is true or false
 * \return _SELF_ Chrom without the elements that overlap the interval
 */
template <class _SELF_,class _BASE_>
long int uGenericNGSChrom<_SELF_,_BASE_>::getDistinctCount(double p_start, double p_end, OverlapType overlap)const
{
    try
    {


        int overlapCount = this->count();
        auto posIter = this->findPrecedingSite(p_start);

        /**<  If no tag leftwise, we start at beginning*/
        if (posIter== this->end())
            posIter= this->begin();
        for (; posIter !=  this->end(); posIter++)
        {
            if (sortGetStart(&(*posIter))> p_end)
                break;

            if (utility::isOverlap(sortGetStart(&(*posIter)), sortGetEnd(&(*posIter)),p_start, p_end,overlap))
            {
                overlapCount--;
            }
        }
        return overlapCount;
    }
    catch (unsorted_throw & e)
    {
#ifdef DEBUG
        std::cerr << "getDistinctCount called on unsorted vector" <<std::endl;
#endif
        throw e;
    }
    catch (ugene_exception_base & e)
    {
#ifdef DEBUG
        std::cerr << "Calling getDistinctCount and you did not provide an appropriate get function" <<std::endl;
#endif
        throw e;
    }
}



/** \brief Split each item into smaller equal size members and replace our vector of items with the new one.
 *
 *      This function will call divideIntoNBin() on each element of the collection, and replace the current elements
 *      with the result of this call. This will invalid most secondary parameters, such as signals, scores and such.
 *
 * \param N int : Number of bins to make
 * \param type SplitType : How to manage splits that are not a multiple of N. Possible : STRICT - IGNORE - FILL - ADD
 * \sa divideIntoNBin
 * \return void
 *
 */
template <class _SELF_,class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::divideItemsIntoNBins(int N, SplitType type)
{
    std::vector<_BASE_> newVector;
    std::vector<_BASE_> tempVector;
    try
    {
        for( _BASE_& x : VecSites)
        {
            tempVector = move(x.divideIntoNBin(N,type));
            newVector.insert( newVector.end(), tempVector.begin(), tempVector.end() );
        }

        VecSites=move(newVector);
    }
    catch(std::exception & e)
    {
#ifdef DEBUG
        std::cerr << "Failed in dividerItemsIntoNBins"<<std::endl;
#endif
        throw e;
    }
}
/** \brief Divide each member into smaller members of size N, sorting the leftover according to our SplitType
 *
 *  This will call divideIntoBinofSize on every element and replace the current elements with the resulting elements.
 *
 *
 * \param N int : Size of the members to build in BP
 * \param type SplitType : How to manage splits that are not a multiple of N. Possible : STRICT - IGNORE - FILL - ADD
 * \sa divideIntoBinofSize
 * \return void
 *
 */
template <class _SELF_,class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::divideItemsIntoBinofSize(int N, SplitType type)
{
    std::vector<_BASE_> newVector;
    std::vector<_BASE_> tempVector;
    try
    {
        for( _BASE_& x : VecSites)
        {
            tempVector =(move(x.divideIntoBinofSize(N,type)));
            for ( unsigned int i = 0; i < tempVector.size(); ++i )
            {
                newVector.push_back(static_cast<_BASE_>( tempVector.at(i) ));
            }
        }
        VecSites=move(newVector);
    }
    catch(std::exception &e)
    {
#ifdef DEBUG
        std::cerr << "Failed in divideItemsIntoBinofSize"<<std::endl;
#endif
        throw e;
    }
}
template <class _SELF_,class _BASE_>
/** \brief Will set scaffold/chr size to be equal to the largest end position of all elements in the collection.
 *  \sa setChromSize
 */
void uGenericNGSChrom<_SELF_,_BASE_>:: inferChrSize()
{
    if (VecSites.size()!=0)
        this->setChromSize(this->maxSite(comparePos)->getEnd());
    else
        this->setChromSize(0);
}

/**<  Wrappers around STL algorithms*/
/** \brief Compute a value for all sites in the chromosome and return the resulting collection
  *
  * This function take a pointer to a functor to perform on all the
  * sites in the collection. Specifically, this will make an std::vector copy
  * of all elements. It will then run the functor using for_each and return the resulting results.
  *
  * \param unary_op UnaryOperation : Unary operation to perform on all the sites of the chromosome
  * \return A collection of values computed on each site by unary_op
  */
template <class _SELF_,class _BASE_>
template<class UnaryOperation>
std::vector<_BASE_> uGenericNGSChrom<_SELF_,_BASE_>::applyAndGetVecData(UnaryOperation unary_op)
{
    std::vector<_BASE_> copyVec(VecSites);
    std::for_each(std::begin(copyVec), std::end(copyVec), unary_op);
    return copyVec;
}

/**<  Wrappers around STL algorithms*/
/** \brief Compute a value for all sites in the chromosome and return the resulting collection
  *
  * This function take a pointer to a functor to perform on all the
  * sites in the collection. Specifically, this will make an std::vector copy
  * of all elements. It will then run the functor using for_each and return the resulting results.
  *
  * \param unary_op UnaryOperation : Unary operation to perform on all the sites of the chromosome
  * \return A collection of values computed on each site by unary_op
  */
template <class _SELF_,class _BASE_>
template<class UnaryOperation>
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::applyAndGetChrom(UnaryOperation unary_op)
{

    _SELF_ returnChrom = this->getCopy();
    returnChrom.applyOnAllSites(unary_op);
    return returnChrom;
}

/** \brief Create a copy of the sites vector, transform it and return the copy
    *
    * This function take a pointer to a function to transform the copied sites
    * vector; It will operator transform() on the elements and return a vector of results.
    * The vector will contain the elements specified by the functor, so may be any value.
    * Please note, the standard requires the functor given to transform() to be side-effect free.
    *
    * The function  passed must have a return value.
    *
    * \exception std::exception : Any exception throw by reserve() or transform()
    * \param unary_op UnaryOperation : Unary operation to perform on the copied sites vector
    * \return A vector of the same type and length as the sites vector but with its sites transformed by unary_op
    */
template <class _SELF_,class _BASE_>
template<class UnaryOperation>
auto uGenericNGSChrom<_SELF_,_BASE_>::computeOnAllSites(UnaryOperation unary_op) -> std::vector<decltype(unary_op(_BASE_()))>
{
    try{
        std::vector<decltype(unary_op(_BASE_()))> results;
        results.reserve(VecSites.size());
        transform(std::begin(VecSites), std::end(VecSites), std::back_inserter(results), unary_op);
        return results;
    }
    catch(...)
    {
        throw;
    }
}
/** \brief Get the sites for which a certain predicate is true
  *
  * This function take a pointer to a predicate function. It return a vector containing
  * the elements of our collection that evalue true when the predicate is applied.
  *
  * The function passed must return a boolean; true if the predicate is true, false otherwise.
  *
  * \param p UnaryPredicate : Unary predicate to evaluate on all sites
  * \return A collection containing all the sites for which the predicate is true
  */
//TODO return chrom?,
template <class _SELF_,class _BASE_>
template<class UnaryPredicate>
std::vector<_BASE_> uGenericNGSChrom<_SELF_,_BASE_>::getSpecificSites(UnaryPredicate pred) const
{
    try
    {
        std::vector<_BASE_> copyVec {};
        // copyVec.reserve(VecSites.size());
        copy_if(std::begin(VecSites), std::end(VecSites), std::back_inserter(copyVec), pred);
        return copyVec;
    }
    catch(std::exception & e)
    {
#ifdef DEBUG
        std::cerr << "Throwing in getSpecificSites()" <<std::endl;
#endif
        throw ugene_exception_base()<<string_error("Caugh std::exception, throwing in getSpecificSites()");
    }
}
/** \brief Remove sites for which the predicate is true.
 *
 *  This function calls remove_if using the supplied predicated, followed by erase(). This will
 *  remove every element that evaluates true from the container  . This will preserve the relative order of the non-removed elements.
 *
 *  The function must return a boolean; true if the predicate is true, false otherwise.
 *
 * \exception std::exception : Any exception throw by erase() or remove_if()
 * \param pred UnaryPredicate Predicate to test, follows standard pattern
 * \return void
 *
 */
template <class _SELF_,class _BASE_>
template<class UnaryPredicate>
void uGenericNGSChrom<_SELF_,_BASE_>::removeSpecificSites(UnaryPredicate pred)
{
    try
    {
        VecSites.erase(std::remove_if(VecSites.begin(), VecSites.end(), pred), VecSites.end());
    }
    catch(...)
    {
#ifdef DEBUG
        std::cerr << "Throwing in removeSpecificSites()" <<std::endl;
#endif
        throw;
    }
}

/** \brief applyOnAllSites: Transform the sites collection by applying a certain function to all sites
    *
    * In-place complement to applyAndGetVecData. This function will run the given functor
    * on every element of the collection, using for_each.
    *
    * If collection is empty, nothing will be done.
    *
    * The function must return void (any other return value will be ignored).
    *
    * \param unary_op UnaryOperation : Unary operation to perform on the sites collection
    * \return unary_op, the operation that was performed on all sites
    */
template <class _SELF_,class _BASE_>
template<class UnaryOperation>
UnaryOperation uGenericNGSChrom<_SELF_,_BASE_>::applyOnAllSites(UnaryOperation f)
{
    try
    {
        if (VecSites.size()>0)
            return for_each(std::begin(VecSites), std::end(VecSites), f);
        else
            return f;
    }
    catch(...)
    {
        throw ;
    }
}

/** \brief applyOnAllSitesConst: Pool the sites collection by applying a certain function to all sites
  *
  *  Const version, see non-const for doc.
  *
  * \sa applyOnAllSites
  */
template <class _SELF_,class _BASE_>
template<class UnaryOperation>
UnaryOperation uGenericNGSChrom<_SELF_,_BASE_>::applyOnAllSites(const UnaryOperation f) const
{
    try
    {
        if (VecSites.size()>0)
            return for_each(std::begin(VecSites), std::end(VecSites), f);
        else
            return f;
    }
    catch(...)
    {
        throw ;
    }
}

/** \brief Accumulate information by querying all sites
  *
  *  Runs accumulate() on the elements of the collection with the given functor. This allows the
  *  querying of every site in a way that returns a single value. ex : adding every elem lenght to
  *  obtain the total lenght of all contigs in the collection.
  *
  *  The function must return the new value of the accumulator.
  *
  * \param binary_op BinaryOperation : Querying function to perform on the sites collection
  * \param init InitialValue The initial value of the "accumulator". Typically 0 if working with an int.
  * \return The information accumulated by querying all the sites
  */
template <class _SELF_,class _BASE_>
template<class BinaryOperation, class InitialValue>
InitialValue uGenericNGSChrom<_SELF_,_BASE_>::accumulateSitesInfo(BinaryOperation binary_op, InitialValue init) const
{
    // Force using sequential version for accumulate as parallel version
    // doesn't work if actual data type of InitialValue and _BASE_ cannot be
    // converted back and forth.
    return __gnu_parallel::accumulate(std::begin(VecSites), std::end(VecSites), init, binary_op, __gnu_parallel::sequential_tag());

}


/** \brief Compute the number of sites for which a certain predicate is true
  *
  * Thie function predicate passed is used to count how many elements in the collection correspond to
  * a given collection.
  *
  * The function must return a boolean; true is the predicate is true, false otherwise.
  *
  * \param p UnaryPredicate : Unary predicate to evaluate on all sites
  * \return The number of sites for which a certain predicate is true
  */
template <class _SELF_,class _BASE_>
template <class UnaryPredicate>
typename std::iterator_traits<typename std::vector<_BASE_>::iterator>::difference_type uGenericNGSChrom<_SELF_,_BASE_>::countSitesWithProperty(UnaryPredicate p) const
{
    try
    {
        return count_if(std::begin(VecSites), std::end(VecSites), p);
    }
    catch(...)
    {
        throw;
    }
}

/**< End wrappers */

/** \brief Sort the sites vector by applying a certain comparison
     *
     * Sort the elements of the collection according to the given binary comparison operator.
     *
     * Additionally, one may provide a pointer to related getters. This enables the use of getSubset()
     * and removeSubset() on the appropriate type of sort. If only get_start_funct is provided, getEnd_funct is set to get_start_funct
     *
     *\
     * \param getStart_funct : function object taking a _BASE_ as member object and returning a value used to sort
     * \param getEnd_funct: function object taking a _BASE_ as member object and returning a value used to break sorting ties.
     * \param comp Compare : Binary comparison operation to perform on the sites collection
     * \return void
     */
template <class _SELF_,class _BASE_>
template <class Compare>
void uGenericNGSChrom<_SELF_,_BASE_>::sortSites(Compare comp,std::function<float(const _BASE_*)> getStart_funct,std::function<float(const _BASE_*)> getEnd_funct)
{
    try
    {
        this->m_isSorted=true;
        sortGetStart=getStart_funct;
        if (getEnd_funct==nullptr)
            sortGetEnd=getStart_funct;
        else
            sortGetEnd= getEnd_funct;
        m_comptFunc=comp;

        return std::sort(std::begin(VecSites), std::end(VecSites), comp);
    }
    catch(...)
    {
        throw;
    }
}

/** \brief Default sort using the start position as a the comparison point
  *
  * \return void
  */
template <class _SELF_,class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::sortSites()
{
    try
    {
        return sortSites(compareStart,&_BASE_::getStart,&_BASE_::getEnd);
    }
    catch (...)
    {
        throw;
    }

}
/** \brief Indicates if the sites collection is sorted according to a certain comparison
  *
  * The function passed must return a boolean: true if the first element is "lower" than the second, false
  * otherwise.
  *
  * \param comp Compare : Binary comparison operation to perform on the sites collection
  * \return true if the sites are sorted, false otherwise.
  */
template <class _SELF_,class _BASE_>
template<class Compare>
bool uGenericNGSChrom<_SELF_,_BASE_>::isSorted(Compare comp) const
{
    return is_sorted(std::begin(VecSites), std::end(VecSites), comp);
}

/** \brief Indicates if the sites collection is sorted ascendingly according the function
  * that was used for the sort (by default: compareStart, which is a sort based on the
  * position of the element on the chromosome).
  *
  * \return true if the sites are sorted, false otherwise.
  */
template <class _SELF_,class _BASE_>
bool uGenericNGSChrom<_SELF_,_BASE_>::isSorted() const
{
    return isSorted(m_comptFunc);
}
/** \brief Find the minimum site according to a certain comparison
  *
  *  The function must return a boolean: true if the first element is "lower" than the second, false otherwise.
  *
  * \param comp Compare : Binary comparison operation to perform on the sites collection
  * \return An iterator to the maximal site
  * \sa maxSite
  * \sa minAndMaxSites
  */
template <class _SELF_,class _BASE_>
template<class Compare>
typename std::vector<_BASE_>::const_iterator uGenericNGSChrom<_SELF_,_BASE_>::minSite(Compare comp) const
{
    try
    {
        return min_element(std::begin(VecSites), std::end(VecSites), comp);
    }
    catch(...)
    {
        throw;
    }
}
/** \brief Find the maximal site according to a certain comparison
  *
  *  The function must return a boolean: true if the first element is "lower" than the second, false otherwise.
  *
  * \param comp Compare : Binary comparison operation to perform on the sites collection
  * \return An iterator to the maximal site
  * \sa minSite
  * \sa minAndMaxSites
  */
template <class _SELF_,class _BASE_>
template<class Compare>
typename std::vector<_BASE_>::const_iterator uGenericNGSChrom<_SELF_,_BASE_>::maxSite(Compare comp) const
{
    try
    {
        return max_element(std::begin(VecSites), std::end(VecSites), comp);
    }
    catch(std::exception & e)
    {
        throw e;
    }
}

/** \brief Find the minimal and maximal sites according to a certain comparison
  *
  * This function take a pointer to a function to find the minimal and
  * maximal sites; this function pointer can either be a) the name of a
  * function taking two sites as parameters, b) a lambda function taking two
  * sites as parameters or c) a member method of a site taking another site
  * as parameter using "mem_fun_ref". In all cases, the function must return
  * a boolean: true if the first element is "lower" than the second, false
  * otherwise.
  *
  * \param comp Compare : Binary comparison operation to perform on the sites collection
  * \return A pair of iterators: the first indicates the minimal site and the second, the maximal site
  */
template <class _SELF_,class _BASE_>
template<class Compare>
std::pair<typename std::vector<_BASE_>::const_iterator, typename std::vector<_BASE_>::const_iterator> uGenericNGSChrom<_SELF_,_BASE_>::minAndMaxSites(Compare comp) const
{
    try
    {
        return minmax_element(std::begin(VecSites), std::end(VecSites), comp);
    }
    catch(std::exception & e)
    {
        throw e;
    }

}







    /** \brief Return a copy of the functor  used to access to current Start value
      *
      *  This will return a copy of the assigned functor to access the start value. Will equal nullptr if not set
      *
      *
      * \return std::function<float(const _BASE_*)> const Copy of the fucntor to access start.
      */
    template <class _SELF_,class _BASE_>
    std::function<float(const _BASE_*)> uGenericNGSChrom<_SELF_,_BASE_>::getStartFunct() const
    {
        return sortGetStart;
    }
    /** \brief Return a copy of the functor  used to access to current end value
    *
    *  This will return a copy of the assigned functor to access the end value. Will equal nullptr if not set
    *
    * \return std::function<float(const _BASE_*)> const Copy of the fucntor to access end.
    */
    template <class _SELF_,class _BASE_>
    std::function<float(const _BASE_*)> uGenericNGSChrom<_SELF_,_BASE_>::getEndFunct() const
    {
        return sortGetEnd;
    }

    /** \brief Return a copy of the comparison functor currently used
     *
     *  Will return a copy of the comparisong functor used for sorting. The comparison functor takes
     *  two _BASE_ elements and returns true or false comparison. Set to nullptr by default
     *
     *
     * \return std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> Comparison functor
     *
     */
    template <class _SELF_,class _BASE_>
    std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> uGenericNGSChrom<_SELF_,_BASE_>::getCompFunct() const
    {
        return m_comptFunc;
    }

    /** \brief Return sorted status of the elements
     *
     * \return bool. True if sorted.
     *
     */
    template <class _SELF_,class _BASE_>
    bool uGenericNGSChrom<_SELF_,_BASE_>::getSortedStatus() const
    {
        return m_isSorted;
    }

    /** \brief  Return number of elements
     *
     * \return int. Number of elements contained
     *
     */
    template <class _SELF_,class _BASE_>
    int uGenericNGSChrom<_SELF_,_BASE_>::count() const
    {
        return VecSites.size();
    }

    /** \brief Return size of scaffold/chrom. 0 if not set
     *
     * \return long long int. Size of the scaffold/chrom.
     *
     */
    template <class _SELF_,class _BASE_>
    long long int uGenericNGSChrom<_SELF_,_BASE_>::getChromSize() const
    {
        return chromSize;
    }

    /** \brief Set scaffold/chrom size
     *
     *  Set the scaffold/chrom size. Must be above 0
     * \param chromS long int. Size to set the scaffold/chrom size to
     * \exception param_throw. Parameter chromS is < then 0.
     * \return void
     *
     */
    template <class _SELF_,class _BASE_>
    void uGenericNGSChrom<_SELF_,_BASE_>::setChromSize(long long int chromS)
    {
        if (chromS<0)
            throw param_throw()<<string_error("failling in setChromSize, value "+utility::to_string(chromS)+" is below 0\n");
        chromSize= chromS;
    }

    /** \brief return name of the scaffold/chrom
     *
     *  Returns a std::string containing the name of the scaffold.
     *  An empty string ("") is a valid name.
     *
     * \return std::string. Name of the scaffold/chrom
     *
     */
    template <class _SELF_,class _BASE_>
    std::string uGenericNGSChrom<_SELF_,_BASE_>::getChr() const
    {
        return chr;
    }

    /** \brief Set the name of the scaffold/chr.
     *
     *   Set the name of the scaffold/chr. This function is included for
     *   working directly with a chrom structure. If working through a Exp structure, this should be set
     *   via a call to the corresponding Exp function, otherwise mapping may be throw off.
     *
     *  \param pChr std::string. Scaffold name to set to.
     *  \return void
     *
     */
    template <class _SELF_,class _BASE_>
    void uGenericNGSChrom<_SELF_,_BASE_>::setChr(std::string pChr)
    {
        chr = move(pChr);
    }


    /** \brief  Return copy of the element at .begin()+position count from iterator
     *
     * \param position int. Position of the _BASE_ in Vecsites to returné (Unrelated to the start and end position of the element)
     * \exception param_throw. Throw if the requested element is an invalid position.
     * \return _BASE_ Copy of the element at the asked for position.
     *
     */
    template <class _SELF_,class _BASE_>
    _BASE_ uGenericNGSChrom<_SELF_,_BASE_>::getSite(long long pPos) const
    {
        try
        {
            return VecSites.at(pPos);
        }
        catch (std::exception &e)
        {
            throw param_throw()<<string_error("Failling in getSite(), out_of_range error");
        }
    }

    /** \brief  Return a vector containing all elements of chrom structure.
     *
     * \return std::vector<_BASE_> std::vector containing all _BASE_ elements in chrom structure
     *
     */
    template <class _SELF_,class _BASE_>
    std::vector<_BASE_> uGenericNGSChrom<_SELF_,_BASE_>::returnVecData() const
    {
        return VecSites;
    }


} // End of namespace NGS
#endif // UFORMATCHROM_H_INCLUDED
