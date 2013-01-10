#ifndef UFORMATCHROM_H_INCLUDED
#define UFORMATCHROM_H_INCLUDED
#include <climits>
#include <iostream>
#include <map>
#include "utility/utility.h"
#include <algorithm>
#include <parallel/numeric>
#include <functional>
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
    bool m_isSorted=true; /*!< If we are in a sorted state or not */
    std::function<float(const _BASE_*)> sortGetStart=nullptr;  /*!<Pointer to the starting sort value */
    std::function<float(const _BASE_*)> sortGetEnd=nullptr ;  /*!< Pointer to the end starting sort value. Typically set to equal Start */
    std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> m_comptFunc=compareStart; /*!< Pointer to sorting function */
    long long int chromSize=0; /*!< Size of the scaffold */

private :
    /**< removeSites overloads */
    void removeSite(int position);
    void removeSite(int start,int end);
    void removeSite(VecGenConstIter position);
    void removeSite(VecGenConstIter start,VecGenConstIter end);


    //TODO move to utility
    template <class Container>
    typename Container::iterator to_mutable_iterator(Container& c, typename Container::const_iterator it)
    {
        return c.begin() + (it - c.begin());
    }

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

    std::vector<long long> returnSiteSizes() const;
    void addDataNoCheck( _BASE_ newSite);
    template<class L, typename S, typename R> friend class uGenericNGSExperiment;

public:

    /**< Constructors */
    uGenericNGSChrom() {};
    uGenericNGSChrom(const std::string & consString):chr(consString) {};
    uGenericNGSChrom(const std::string & consString, long int size);
    uGenericNGSChrom(const std::vector<_BASE_> &);

    uGenericNGSChrom& operator=(const uGenericNGSChrom& copFrom);
    uGenericNGSChrom(const uGenericNGSChrom&);

    /**< Declared functions */
    unsigned long long avgSiteSize() const;
    unsigned long long minSiteSize() const;
    unsigned long long maxSiteSize() const;
    unsigned long long sumSiteSize() const;
    void printStats(std::ostream& out) const;

    /** \brief Default virtual destructor of uGenericNGSChrom
     */
    virtual ~uGenericNGSChrom<_SELF_,_BASE_> ()
    {
        ;
    }

    void inferChrSize();

    //TODO MAKE THIS WORK
    long int countUnique() const;

    /**< In place */
    void divideItemsIntoNBins(int N, SplitType type=SplitType::STRICT);
    void divideItemsIntoBinofSize(int N, SplitType type=SplitType::STRICT);

    /**< Find according to sort value */
    VecGenConstIter findPrecedingSite(const float position) const;
    VecGenConstIter findNextSite(const float position) const;
    /**< Functions to create and add items to our chrom */
    template <class _OTHER_>
    _BASE_ generateRandomSite(const int size, std::mt19937& engine, const _OTHER_ &exclList, const int sigma=0, const std::string ID="") const;
    _BASE_ generateRandomSite(const int size, std::mt19937& engine, const int sigma=0, const std::string ID="") const;

    template <class _OTHER_>
    void addNRandomSite(const int size, const int n, std::mt19937& engine, const _OTHER_ &exclList, const int sigma=0, const std::string ID="");
    void addNRandomSite(const int size, const int n, std::mt19937& engine, const int sigma=0, const std::string ID="");

    template <class _OTHER_>
    _SELF_ getOverlapping(_OTHER_ &compareChr,OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;
    template <class _OTHER_>
    _SELF_ getNotOverlapping(_OTHER_ &compareChr,OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;

    /**<  */
    _SELF_ getDistinct(float p_start, float p_end, OverlapType options=OverlapType::OVERLAP_PARTIAL) const;

    /**< Functions to manipulate generically ranges of our elements */
    _SELF_ getSubset(float p_start, float p_end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;
    _SELF_ removeSubset(float p_start, float p_end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);

    long int getSubsetCount(float p_start, float p_end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL)const;

    void addData(const _BASE_ & newSite);

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

    /** \brief Return a copy of the functor  used to access to current Start value
     *
     *  This will return a copy of the assigned functor to access the start value. Will equal nullptr if not set
     *
     *
     * \return std::function<float(const _BASE_*)> const Copy of the fucntor to access start.
     */
    std::function<float(const _BASE_*)> getStartFunct() const
    {
        return sortGetStart;
    }
    /** \brief Return a copy of the functor  used to access to current end value
    *
    *  This will return a copy of the assigned functor to access the end value. Will equal nullptr if not set
    *
    * \return std::function<float(const _BASE_*)> const Copy of the fucntor to access end.
    */
    std::function<float(const _BASE_*)> getEndFunct() const
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
    std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> getCompFunct() const
    {
        return m_comptFunc;
    }

    /** \brief Return sorted status of the elements
     *
     * \return bool. True if sorted.
     *
     */
    bool getSortedStatus() const
    {
        return m_isSorted;
    }

    /** \brief  Return number of elements
     *
     * \return int. Number of elements contained
     *
     */
    int count() const
    {
        return VecSites.size();
    };

    /** \brief Return size of scaffold/chrom. 0 if not set
     *
     * \return long long int. Size of the scaffold/chrom.
     *
     */
    long long int getChromSize() const
    {
        return chromSize;
    };

    /** \brief Set scaffold/chrom size
     *
     *  Set the scaffold/chrom size. Must be above 0
     * \param chromS long int. Size to set the scaffold/chrom size to
     * \exception param_throw. Parameter chromS is < then 0.
     * \return void
     *
     */
    void setChromSize(long long int chromS)
    {
        {
            if (chromS<0)
                throw param_throw()<<string_error("failling in setChromSize, value "+utility::to_string(chromS)+" is below 0\n");
            chromSize= chromS;
        }
    };

    /** \brief return name of the scaffold/chrom
     *
     *  Returns a std::string containing the name of the scaffold.
     *  An empty string ("") is a valid name.
     *
     * \return std::string. Name of the scaffold/chrom
     *
     */
    std::string getChr() const
    {
        return chr;
    };

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
    void setChr(std::string pChr)
    {
        chr = move(pChr);
    };


    /** \brief  Return copy of the element at .begin()+position count from iterator
     *
     * \param position int. Position of the _BASE_ in Vecsites to returné (Unrelated to the start and end position of the element)
     * \exception param_throw. Throw if the requested element is an invalid position.
     * \return _BASE_ Copy of the element at the asked for position.
     *
     */
    _BASE_ getSite(int position)
    {
        try
        {
            return VecSites.at(position);
        }
        catch (std::exception &e)
        {
            throw param_throw()<<string_error("Failling in getSite(), out_of_range error");
        }
    };

    /** \brief  Return a vector containing all elements of chrom structure.
     *
     * \return std::vector<_BASE_> std::vector containing all _BASE_ elements in chrom structure
     *
     */
    std::vector<_BASE_> returnVecData() const
    {
        return VecSites;
    };

    /**<  Wrappers around STL algorithms*/
    /** \brief Compute a value for all sites in the chromosome and return the resulting collection
      *
      * This function take a pointer to a functor to perform on all the
      * sites in the collection. Specifically, this will make an std::vector copy
      * of all elements. It will then run the functor using for_each and return the resulting results.
      *
      *
      * The function pointer can either be a)
      * the name of a function taking a site by reference, b) a lambda
      * function taking a site by reference or c) a member method of a
      * site using "mem_fun_ref" or "ref". In all cases, the function must return a
      * non void value.
      *
      * \param unary_op UnaryOperation : Unary operation to perform on all the sites of the chromosome
      * \return A collection of values computed on each site by unary_op
      */
    template<class UnaryOperation>
    std::vector<_BASE_> applyAndGetVecData(UnaryOperation unary_op)
    {
        std::vector<_BASE_> copyVec(VecSites);
        for_each(begin(copyVec), end(copyVec), unary_op);
        return copyVec;
    }

    /** \brief Create a copy of the sites vector, transform it and return the copy
      *
      * This function take a pointer to a function to transform the copied sites
      * vector; It will operator transform() on the elements and return a vector of results.
      * The vector will contain the elements specified by the functor, so may be any value.
      * Please note, the standard requires the functor given to transform() to be side-effect free.
      *
      *
      * The functor can either be a) the name of a function
      * taking a site by reference, b) a lambda function taking a site by
      * reference or c) a member method of a site using "mem_fun_ref" or "ref".
      * The function must have a return value.
      *
      * \exception std::exception : Any exception throw by reserve() or transform()
      * \param unary_op UnaryOperation : Unary operation to perform on the copied sites vector
      * \return A vector of the same type and length as the sites vector but with its sites transformed by unary_op
      */
    //TODO send back a chrom?
    template<class UnaryOperation>
    auto computeOnAllSites(UnaryOperation unary_op) -> std::vector<decltype(unary_op(_BASE_()))>
    {
        try{
            std::vector<decltype(unary_op(_BASE_()))> results;
            results.reserve(VecSites.size());
            transform(std::begin(VecSites), std::end(VecSites), std::back_inserter(results), unary_op);
            return results;
        }
        catch(std::exception &e)
        {
            throw e;
        }
    }

    /** \brief Get the sites for which a certain predicate is true
      *
      * This function take a pointer to a predicate function. It return a vector containing
      * the elements of our collection that evalue true when the predicate is applied.
      *
      *
      * The functor can either be a) the name of a function taking a site as
      * parameter, b) a lambda function taking a site as parameter or c) a
      * member method of a site using "mem_fun_ref" or "ref". In all cases, the function
      * must return a boolean; true if the predicate is true, false otherwise.
      *
      * \param p UnaryPredicate : Unary predicate to evaluate on all sites
      * \return A collection containing all the sites for which the predicate is true
      */
    //TODO return chrom?,
    template<class UnaryPredicate>
    std::vector<_BASE_> getSpecificSites(UnaryPredicate pred) const
    {
        try
        {
            std::vector<_BASE_> copyVec {};
            copyVec.reserve(VecSites.size());
            copy_if(std::begin(VecSites), std::end(VecSites), std::back_inserter(copyVec), pred);
            return copyVec;
        }
        catch(std::exception & e)
        {
#ifdef DEBUG
            std::cerr << "Throwing in getSpecificSites()" <<std::endl;
#endif
            throw e;
        }
    }
    template<class UnaryPredicate>
    /** \brief Remove sites for which the predicate is true.
     *
     *  This function calls remove_if using the supplied predicated, followed by erase(). This will
     *  remove every element that evaluates true from the container  . This will preserve the relative order of the non-removed elements.
     *
     *  The functor can either be a) the name of a function taking a site as
     *  parameter, b) a lambda function taking a site as parameter or c) a
     *  member method of a site using "mem_fun_ref" or "ref". In all cases, the function
     *  must return a boolean; true if the predicate is true, false otherwise.
     *
     * \exception std::exception : Any exception throw by erase() or remove_if()
     * \param pred UnaryPredicate Predicate to test, follows standard pattern
     * \return void
     *
     */
    void removeSpecificSites(UnaryPredicate pred)
    {
        try
        {
            VecSites.erase(std::remove_if(VecSites.begin(), VecSites.end(), pred), VecSites.end());
        }
        catch(std::exception & e)
        {
#ifdef DEBUG
            std::cerr << "Throwing in removeSpecificSites()" <<std::endl;
#endif
            throw e;
        }
    }

    /** \brief applyOnAllSites: Transform the sites collection by applying a certain function to all sites
      *
      * In-place complement to applyAndGetVecData. This function will run the given functor
      * on every element of the collection, using for_each.
      *
      * If collection is empty, nothing will be done.
      *
      * The functor pointer can either be a) the name of a function
      * taking a site by reference, b) a lambda function taking a site by
      * reference or c) a member method of a site using "mem_fun_ref". In all
      * cases, the function must return void (any other return value will be
      * ignored).
      *
      * \param unary_op UnaryOperation : Unary operation to perform on the sites collection
      * \return unary_op, the operation that was performed on all sites
      */
    template<class UnaryOperation>
    UnaryOperation applyOnAllSites(UnaryOperation f)
    {
        try
        {
            if (VecSites.size()>0)
                return for_each(std::begin(VecSites), std::end(VecSites), f);
            else
                return f;
        }
        catch(std::exception & e)
        {
            throw e;
        }
    }


    /** \brief applyOnAllSitesConst: Pool the sites collection by applying a certain function to all sites
      *
      *  Const version, see non-const for doc.
      *
      * \sa applyOnAllSites
      */
    template<class UnaryOperation>
    UnaryOperation applyOnAllSites(const UnaryOperation f) const
    {
        try
        {
            return for_each(std::begin(VecSites), std::end(VecSites), f);
        }
        catch(std::exception & e)
        {
            throw e;
        }
    }

    /** \brief Accumulate information by querying all sites
      *
      *  Runs accumulate() on the elements of the collection with the given functor. This allows the
      *  querying of every site in a way that returns a single value. ex : adding every elem lenght to
      *  obtain the total lenght of all contigs in the collection.
      *
      *
      * The function pointer can either be a) the name of a function taking two
      * parameters, an accumulator and a site or b) a lambda function taking two
      * parameters, an accumulator and a site. In all cases, the function must
      * return the new value of the accumulator.
      *
      * \param binary_op BinaryOperation : Querying function to perform on the sites collection
      * \param init InitialValue The initial value of the "accumulator". Typically 0 if working with an int.
      * \return The information accumulated by querying all the sites
      */
    template<class BinaryOperation, class InitialValue>
    InitialValue accumulateSitesInfo(BinaryOperation binary_op, InitialValue init) const
    {
        // Force using sequential version for accumulate as parallel version
        // doesn't work if actual data type of InitialValue and _BASE_ cannot be
        // converted back and forth.
        return __gnu_parallel::accumulate(std::begin(VecSites), std::end(VecSites), init, binary_op, __gnu_parallel::sequential_tag());
    }
    /**< End STL wrappers */


    /** \brief Sort the sites vector by applying a certain comparison
      *
      * Sort the elements of the collection according to the given binary comparison operator.
      *
      * Addtionally, one may provide a pointer to related getters. This enables the use of getSubset()
      * and removeSubset() on the appropriate type of sort. If only get_start_funct is provided, getEnd_funct is set to get_start_funct
      *
      *\
      * \param getStart_funct : function object taking a _BASE_ as member object and returning a value used to sort
      * \param getEnd_funct: function object taking a _BASE_ as member object and returning a value used to break sorting ties.
      * \param comp Compare : Binary comparison operation to perform on the sites collection
      * \return void
      */

    template<class Compare>
    void sortSites(Compare comp,std::function<float(const _BASE_*)> getStart_funct=nullptr,std::function<float(const _BASE_*)> getEnd_funct=nullptr)
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
        catch(std::exception &e)
        {
            throw e;
        }
    }

    /** \brief Default sort using the start position as a the comparison point
      *
      * \return void
      */
    void sortSites()
    {
        try
        {
            return sortSites(compareStart,&_BASE_::getStart,&_BASE_::getEnd);
        }
        catch (std::exception & e )
        {
            throw e;
        }

    }
    /** \brief Indicates if the sites collection is sorted according to a certain comparison
      *
      * This function take a pointer to a function to determine if the the sites
      * collection is sorted; this function pointer can either be a) the name of a
      * function taking two sites as parameters, b) a lambda function taking two
      * sites as parameters or c) a member method of a site taking another site
      * as parameter using "mem_fun_ref". In all cases, the function must return
      * a boolean: true if the first element is "lower" than the second, false
      * otherwise.
      *
      * \param comp Compare : Binary comparison operation to perform on the sites collection
      * \return true if the sites are sorted, false otherwise.
      */
    template<class Compare>
    bool isSorted(Compare comp) const
    {
        return is_sorted(std::begin(VecSites), std::end(VecSites), comp);
    }

    /** \brief Indicates if the sites collection is sorted ascendingly according to
      * their start position
      *
      * \return true if the sites are sorted, false otherwise.
      */
    bool isSorted() const
    {
        return isSorted(compareStart);
    }

    /** \brief Find the minimal site according to a certain comparison
      *
      * This function take a pointer to a function to find the minimal site;
      * this function pointer can either be a) the name of a function taking two
      * sites as parameters, b) a lambda function taking two sites as parameters
      * or c) a member method of a site taking another site as parameter using
      * "mem_fun_ref". In all cases, the function must return a boolean: true if
      * the first element is "lower" than the second, false otherwise.
      *
      * \param comp Compare : Binary comparison operation to perform on the sites collection
      * \return An iterator to the minimal site
      */
    template<class Compare>
    VecGenConstIter minSite(Compare comp) const
    {
        try
        {
            return min_element(std::begin(VecSites), std::end(VecSites), comp);
        }
        catch(std::exception & e)
        {
            throw e;
        }
    }

    /** \brief Find the maximal site according to a certain comparison
      *
      * This function take a pointer to a function to find the maximal site;
      * this function pointer can either be a) the name of a function taking two
      * sites as parameters, b) a lambda function taking two sites as parameters
      * or c) a member method of a site taking another site as parameter using
      * "mem_fun_ref". In all cases, the function must return a boolean: true if
      * the first element is "lower" than the second, false otherwise.
      *
      * \param comp Compare : Binary comparison operation to perform on the sites collection
      * \return An iterator to the maximal site
      */
    template<class Compare>
    VecGenConstIter maxSite(Compare comp) const
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
    template<class Compare>
    std::pair<VecGenConstIter, VecGenConstIter> minAndMaxSites(Compare comp) const
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

    /** \brief Compute the number of sites for which a certain predicate is true
      *
      * This function is used to count how many elements in the collection correspond to
      * a given collection.
      *
      * This function pointer can either be a) the name of a function taking a site as
      * parameter, b) a lambda function taking a site as parameter or c) a
      * member method of a site using "mem_fun_ref" or "ref". In all cases, the function
      * must return a boolean; true is the predicate is true, false otherwise.
      *
      * \param p UnaryPredicate : Unary predicate to evaluate on all sites
      * \return The number of sites for which a certain predicate is true
      */
    template <class UnaryPredicate>
    typename std::iterator_traits<VecGenIter>::difference_type
    countSitesWithProperty(UnaryPredicate p) const
    {
        try
        {
            return count_if(begin(VecSites), end(VecSites), p);
        }
        catch(std::exception & e)
        {
            throw e;
        }
    }
    /**< End STL wrappers */


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
        VecSites.push_back(std::move(newSite));
    }
    catch(ugene_exception_base & e)
    {
        throw e;
    }
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
        m_isSorted=false;
        VecSites.erase(VecSites.begin()+(position));
    }
    catch(std::exception & e)
    {
        throw e;
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
    catch(std::exception & e)
    {
        throw e;
    }
}
template <class _SELF_, class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::removeSite(VecGenConstIter start,VecGenConstIter end)
{
    try
    {
        m_isSorted=false;
        /**< According to C++11 standard, const iterator shoudl be allowed in erase
        However, implementation does not seem to have caught up, hence this patch
         */

        VecSites.erase(to_mutable_iterator(VecSites,start),to_mutable_iterator(VecSites,end));
    }
    catch(std::exception &e )
    {
        throw e;
    }
}
template <class _SELF_, class _BASE_>
void uGenericNGSChrom<_SELF_,_BASE_>::removeSite(VecGenConstIter position)
{
    try
    {
        m_isSorted=false;
        VecSites.erase(to_mutable_iterator(VecSites,position),to_mutable_iterator(VecSites,position));
    }
    catch(std::exception &e)
    {
        throw e;
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
        _SELF_ emptyExcl;
        return generateRandomSite(size,engine,emptyExcl,sigma,ID);
    }
    catch(std::exception &e)
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
(const int size_,std::mt19937& engine,const _OTHER_ &exclList, const int sigma, const std::string ID) const
{
    //TODO Sanity check here to make sure it is possible to generate the asked for tag.
    try
    {
        _BASE_ returnTag;

        bool found=false;
        int size = size_;

        int max = this->getChromSize();

        while (!found)
        {
            int max = this->getChromSize();
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
                    temptag.setEnd(center+shift);
                    temptag.setStart(center-shift);
                    temptag.setChr(this->getChr());
                    //TODO DO MAKE THIS WORK
                    temptag.setName(ID);
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
    catch(std::exception & e)
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
    catch(std::exception &e)
    {
        throw e;
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
    catch(std::exception &e)
    {
        throw e;
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
    try
    {
        if (this->count() == 0)
            return ULONG_MAX;


        return minSite(compareLenght)->getLenght();
    }
    catch(std::exception & e)
    {
        throw e;
    }
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
    long int count=0;
    int current;
    std::map<int, int >  myUniqueMap;

    typename std::vector<_BASE_>::iterator iterVec;

    for (iterVec = VecSites.begin() ; iterVec!= VecSites.end(); iterVec++)

//Output functions

    template <class _SELF_,class _BASE_>
    void uGenericNGSChrom<_SELF_,_BASE_>::printStats(std::ostream& out) const
    {
        current= iterVec.getStart();
        if (myUniqueMap.count(current)==0)
            myUniqueMap.insert( std::pair<int,int>(current,current));
    }
    count=myUniqueMap.size();

    return count;
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
 * \param position const int& value to evaluate from
 * \return typename std::vector<_BASE_>::const_iterator Constant iterator pointing to the item
 *
 */
template <class _SELF_,class _BASE_>
typename std::vector<_BASE_>::const_iterator uGenericNGSChrom<_SELF_,_BASE_>::findPrecedingSite(const float position) const
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
        auto lower = std::lower_bound(VecSites.begin(), VecSites.end(), position, comp);

        /**< If no result, or result is our first item */
        if (lower==VecSites.end()||(lower==VecSites.begin()))
            return VecSites.end();

        /**<Return item precedes and as such is LESS then position  */
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
template <class _SELF_,class _BASE_>
/** \brief Find first element after a given value for our current sort type. Data must be sorted
 *
 *  This function will find the element after the given value as valid to the current sort condition. By default,
 *  the sort condition is based on the start position of every element, but proper use of the sorting function can modify this condition.
 *
 *  If called on an unsorted collection, this function will raise an exception.
 *
 * \exception unsorted_throw : Will throw if the collection is unsorted
 * \exception ugene_exception_base : Will throw if the proper getters where not set
 * \param position const int& value to evaluate from
 * \return typename std::vector<_BASE_>::const_iterator Constant iterator pointing to the value
 *
 */
typename std::vector<_BASE_>::const_iterator uGenericNGSChrom<_SELF_,_BASE_>::findNextSite(const float position) const
{
    try
    {
        /**< If unsorted, fail */

        if (VecSites.size()==0)
            return VecSites.end();

        if ((m_isSorted==false)||(sortGetStart==nullptr)||(sortGetEnd==nullptr))
            throw ugene_exception_base();
        /**< Return true comparitor if item1 smaller then item 2 */
        auto comp = [&] (const float &item1, const _BASE_ &item2)
        {
            return item1< sortGetStart(&item2);
        };

        /**< Compare, sort Value */
        auto upper = std::upper_bound(VecSites.begin(), VecSites.end(), position, comp);

        /**< If no result, or result is our first item */
        if (upper==VecSites.end())
            return VecSites.end();

        /**<Return the item greater then value*/
        return (upper);
    }
    catch (std::exception & e)
    {
#ifdef DEBUG
        std::cerr << "Calling findNextSite on unsorted vector or you did not provide an approriate get function" <<std::endl;
        std::cerr << "is Nullprt Start "<< (sortGetStart==nullptr) <<std::endl;
        std::cerr << "is Nullprt End "<< (sortGetEnd==nullptr) <<std::endl;
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
 * \param start float start interval value
 * \param end float end interval value
 * \param overlap OverlapType OverlapType from Enum
 * \return int Number of elements in our range
 *
 */
template <class _SELF_,class _BASE_>
long int uGenericNGSChrom<_SELF_,_BASE_>::getSubsetCount(float p_start, float p_end, OverlapType overlap) const
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
    catch (ugene_exception_base & e )
    {
        throw e;
    }
    catch(std::exception & e)
    {
        //std::cerr << " Catching in getSubsetCount" <<std::endl;
        throw e;
    }
}


/** \brief Return a subset of our data that overlaps range start/end, based on current sort type.
 *
 * \param start float: Start of interval
 * \param end float: End of interval
 * \param overlap OverlapType: Type of overlap
 * \return uGenericNGSChrom<_BASE_>: Chrom structure containing our element subset
 *
 */
template <class _SELF_,class _BASE_>
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::getSubset(float p_start, float p_end, OverlapType overlap) const
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
    catch (ugene_exception_base & e )
    {
        throw e;
    }
    catch(std::exception & e)
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
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::removeSubset(float p_start, float p_end, OverlapType overlap)
{
    try
    {
        _SELF_ returnChrom;
        returnChrom.setChr(this->getChr());
        std::vector<int> erasePositions;

        auto pos = this->findPrecedingSite(p_start);

        /**<  If no tag leftwise, we start at beginning*/
        if (pos==this->end())
            pos=this->begin();

        auto delPos=pos;
        for (; pos != this->end(); pos++)
        {
            if (sortGetStart(&(*pos))> end)
                break;
            /**< When we find a valid element, go back one step and erase th element */
            if (utility::isOverlap(sortGetStart(&(*pos)), sortGetEnd(&(*pos)),p_start, p_end,overlap))
            {
                returnChrom.addDataNoCheck(*pos);
                pos--;
                delPos=pos;
                this->removeSite( to_mutable_iterator(VecSites,delPos ));
            }
        }

        return returnChrom;
    }
    catch(std::exception & e)
    {
        throw e;
    }

}
/**< Return elements of A that overlap B */
template <class _SELF_,class _BASE_>
template <class _OTHER_>
/** \brief Wrapper function that returns a chrom structure containing the elements of that overlap another chrom structur
 *
 *      This function return a collection. This collection contains every element of THIS that overlaps an element of compareChr. This comparison
 *      is always based on genomic positions ( start end )
 *
 * \param _OTHER_ & compareChr : A compatible chrom collection
 * \param OverlapType overlap  : type of overlap
 * \return Chrom collection containing all overlapping elements.
 * \sa getNotOverlapping
 */

_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::getOverlapping(_OTHER_ &compareChr,OverlapType overlap) const
{
    try
    {
        _SELF_ returnChr;
        for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
        {
            for(auto compit= compareChr.begin(); compit!=compareChr.end(); compit++)
            {
                if (utility::isOverlap(it->getStart(), it->getEnd(),compit->getStart(),compit->getEnd(),overlap))
                {
                    returnChr.addDataNoCheck(*it);
                    break;
                }
            }
        }
        return returnChr;
    }
    catch(std::exception & e)
    {
        throw e;
    }

}


/**< Return elements of A that overlap B */
template <class _SELF_,class _BASE_>
template <class _OTHER_>
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
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::getNotOverlapping(_OTHER_ &compareChr,OverlapType overlap)const
{
    try
    {
        _SELF_ returnChr;
        bool add=true;
        for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
        {
            for(auto compit= compareChr.begin(); compit!=compareChr.end(); compit++)
            {
                if (utility::isOverlap(it->getStart(), it->getEnd(),compit->getStart(),compit->getEnd()))
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
//TODO Test this?
template <class _SELF_,class _BASE_>
/** \brief
 *
 * \param p_start float
 * \param p_end float
 * \param overlap OverlapType
 * \return _SELF_
 *
 */
_SELF_ uGenericNGSChrom<_SELF_,_BASE_>::getDistinct(float p_start, float p_end, OverlapType overlap) const
{
    try
    {
        _SELF_ returnChrom;
        returnChrom= *this;
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
    catch(std::exception & e)
    {

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
//TODO check this again
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
    this->setChromSize(this->maxSite(comparePos)->getEnd());
}

} // End of namespace NGS
#endif // UFORMATCHROM_H_INCLUDED
