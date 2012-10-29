#ifndef UFORMATCHROM_H_INCLUDED
#define UFORMATCHROM_H_INCLUDED
#include <climits>
#include <iostream>
#include <map>
#include "utility/utility.h"
#include <algorithm>
#include <parallel/numeric>
#include <functional>
namespace NGS {
template<typename _BASE_>
class uGenericNGSChrom
{
    static_assert(
        std::is_convertible<_BASE_, uGenericNGS>::value,
        "The type does not inherit from uGenericNGS."
    );
    typedef std::vector<_BASE_> VecGenericNGS;
    typedef typename std::vector<_BASE_>::iterator VecGenIter;
    typedef typename std::vector<_BASE_>::const_iterator VecGenConstIter;

protected:

    std::vector<_BASE_> VecSites{};
    std::string chr="";

private :
    /**< removeSites overloads */
    void removeSite(int position);
    void removeSite(int start,int end);
    void removeSite(VecGenIter position);
    void removeSite(VecGenIter start,VecGenIter end);

    /**< Should not be necesary in C++ 11? */
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

    /**< Pointers to our functions and determines if sorted */
    bool m_isSorted=true;
    std::function<int(const _BASE_*)> sortGetStart=nullptr;
    std::function<int(const _BASE_*)> sortGetEnd=nullptr ;
    std::function<bool(const _BASE_ &item1, const _BASE_ &item2)> m_comptFunc=nullptr;
    long int chromSize=0;



    //TODO Erase these?
    std::vector<long long> returnSiteSizes() const;

    unsigned long long avgSiteSize() const;
    unsigned long long minSiteSize() const;
    unsigned long long maxSiteSize() const;
    unsigned long long sumSiteSize() const;

    template<typename S, typename R> friend class uGenericNGSExperiment;

public:

    void inferChrSize();
    virtual ~uGenericNGSChrom<_BASE_> ()
    {;}

    /**< Write functions */
    void outputBedFormat(std::ostream& out);
    void printStats(std::ostream& out) const;

    /**< In place */
    void divideItemsIntoNBins(int N, SplitType type=SplitType::STRICT);
    void divideItemsIntoBinofSize(int N, SplitType type=SplitType::STRICT);


    /**< Find according to sort value */
    VecGenConstIter findPrecedingSite(const int & position) const;
    VecGenConstIter findNextSite(const int & position) const;
    /**< Functions to create and add items to our chrom */
    template <class T2>
    _BASE_ generateRandomSite(const int size, std::mt19937& engine, const uGenericNGSChrom<T2> &exclList, const int sigma=0, const std::string ID="") const;
    _BASE_ generateRandomSite(const int size, std::mt19937& engine, const int sigma=0, const std::string ID="") const;

    template <class T2>
    void addNRandomSite(const int size, const int n, std::mt19937& engine, const uGenericNGSChrom<T2> &exclList, const int sigma=0, const std::string ID="");
    void addNRandomSite(const int size, const int n, std::mt19937& engine, const int sigma=0, const std::string ID="");


    template <class T2>
    uGenericNGSChrom<_BASE_> getOverlapping(uGenericNGSChrom<T2> &compareExp,OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;
    template <class T2>
    uGenericNGSChrom<_BASE_> getNotOverlapping(uGenericNGSChrom<T2> &compareExp,OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;
    uGenericNGSChrom<_BASE_> getDistinct(std::string chr, int start, int end, OverlapType options=OverlapType::OVERLAP_PARTIAL);

    /**< Functions to manipulate generically ranges of our elements */
    uGenericNGSChrom<_BASE_> getSubset(int start, int end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;
    uGenericNGSChrom<_BASE_> removeSubset(int start, int end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);

    void addSite( _BASE_ newSite);
    int getSubsetCount(int start, int end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL)const;


    bool getSortedStatus() const
    {
        return m_isSorted;
    }

    /**< Return number of elements */
    int count() const
    {
        return VecSites.size();
    };
    /**< Return size */
    int getChromSize() const
    {
        return chromSize;
    };

    /** \brief Set element collection max size
     *
     * \param chromS long int size to set as
     * \return void
     *
     */
    void setChromSize(long int chromS)
    {
        try
        {
            if (chromS<0)
                throw param_throw()<<string_error("failling in setChromSize, value "+std::to_string(chromS)+" is below 0\n");
            chromSize= chromS;
        }
        catch(std::exception e)
        {
            throw e;
        }
    };
    /**< return name of the element. */
    std::string getChr() const
    {
        return chr;
    };

    /**< Set name of the element. Be careful about this, can throw of Experiment mapping */
    void setChr(std::string pchr)
    {
        chr = move(pchr);
    };

    /**< Return copy of the element at .begin()+position count from iterator */
    _BASE_ getSite(int position)
    {
        try
        {
            return VecSites.at(position);
        }
        catch (std::exception)
        {
            throw;
        }
    };

    /**< Return a vector containing all elements. */
    std::vector<_BASE_> returnVecData()
    {
        return VecSites;
    };

    /**<  Wrappers around the STL algorithms*/
    /** \brief Compute a value for all sites in the chromosome and return the resulting collection
      *
      * This function take a pointer to a function to perform on all the
      * sites in the collection; this function pointer can either be a)
      * the name of a function taking a site by reference, b) a lambda
      * function taking a site by reference or c) a member method of a
      * site using "mem_fun_ref". In all cases, the function must return a
      * non void value.
      *
      * \param unary_op UnaryOperation : Unary operation to perform on all the sites of the chromosome
      * \return A collection of values computed on each site by unary_op
      */
//TODO erase
    /*   template<class UnaryOperation>
       std::vector<_BASE_> transformAndGetVecData(UnaryOperation unary_op)
       {
           std::vector<_BASE_> copyVec(VecSites);
           for_each(begin(copyVec), end(copyVec), unary_op);
           return copyVec;
       } */

    /** \brief Create a copy of the sites vector, transform it and return the copy
      *
      * This function take a pointer to a function to transform the copied sites
      * vector; this function pointer can either be a) the name of a function
      * taking a site by reference, b) a lambda function taking a site by
      * reference or c) a member method of a site using "mem_fun_ref". In all
      * cases, the function must return void (any other return value will be
      * ignored).
      *
      * \param unary_op UnaryOperation : Unary operation to perform on the copied sites vector
      * \return A vector of the same type and length as the sites vector but with its sites transformed by unary_op
      */
    //TODO send back a chrom&
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
      * This function take a pointer to a predicate function; this function
      * pointer can either be a) the name of a function taking a site as
      * parameter, b) a lambda function taking a site as parameter or c) a
      * member method of a site using "mem_fun_ref". In all cases, the function
      * must return a boolean; true is the predicate is true, false otherwise.
      *
      * \param p UnaryPredicate : Unary predicate to evaluate on all sites
      * \return A collection containing all the sites for which the predicate is true
      */
    //TODO return chrom,
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


    /** \brief Transform the sites collection by applying a certain function to all sites
      *
      * This function take a pointer to a function to transform the sites
      * collection; this function pointer can either be a) the name of a function
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
        {   if (VecSites.size()>0)
                return for_each(std::begin(VecSites), std::end(VecSites), f);
            else
                return f;
        }
        catch(std::exception & e)
        {
            throw e;
        }
    }

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
          * This function take a pointer to a function to accumulate some information;
          * this function pointer can either be a) the name of a function taking two
          * parameters, an accumulator and a site or b) a lambda function taking two
          * parameters, an accumulator and a site. In all cases, the function must
          * return the new value of the accumulator.
          *
          * \param binary_op BinaryOperation : Querying function to perform on the sites collection
          * \param init InitialValue The initial value of the "accumulator"
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
          * This function take a pointer to a function to sort the sites collection;
          * this function pointer can either be a) the name of a function taking two
          * sites as parameters, b) a lambda function taking two sites as parameters
          * or c) a member method of a site taking another site as parameter using
          * "mem_fun_ref". In all cases, the function must return a boolean: true if
          * the first element is "lower" than the second, false otherwise.
          *
          *
          * Addtionally, one may provide a pointer to related getters. This enables the use of getSubset()
          * and removeSubset() on the appropriate type of sort
          *
          *
          * \param comp Compare : Binary comparison operation to perform on the sites collection
          * \return void
          */

        template<class Compare>
        void sortSites(Compare comp,std::function<int(const _BASE_*)> getStart_funct=nullptr,std::function<int(const _BASE_*)> getEnd_funct=nullptr)
        {
            try
            {
                this->m_isSorted=true;
                sortGetStart=getStart_funct;
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
          * This function take a pointer to a predicate function; this function
          * pointer can either be a) the name of a function taking a site as
          * parameter, b) a lambda function taking a site as parameter or c) a
          * member method of a site using "mem_fun_ref". In all cases, the function
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


        /**< Constructors */
        uGenericNGSChrom() {};
        uGenericNGSChrom(std::string consString):chr(consString) {};
        uGenericNGSChrom(std::string consString, long int size);

        uGenericNGSChrom(std::vector<_BASE_>);


        long int countUnique() const;

        /**< Public iterators */
        auto begin()const->decltype(VecSites.cbegin())
        {
            return VecSites.cbegin();
        };
        auto end()const->decltype(VecSites.cend())
        {
            return VecSites.cend();
        };

          /**< Private iterators */
            auto begin()->decltype(VecSites.begin())
            {
                return VecSites.begin();
            };
            auto end()->decltype(VecSites.end())
            {
                return VecSites.end();
            };

    };



    /** \brief Construct with name and size
     *
     * \param consString std::string
     * \param size long int
     *
     */
    template <class _BASE_>
    uGenericNGSChrom<_BASE_>::uGenericNGSChrom(std::string consString, long int size):chr(consString)
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

    template <class _BASE_>
    /** \brief add a new element to our chrom, throw out_of_mem if impossible
     *
     * \param newSite _BASE_ Ellement to add
     *
     */
    void uGenericNGSChrom<_BASE_>::addSite(_BASE_ newSite)
    {
        try
        {
            m_isSorted=false;
            VecSites.push_back(std::move(newSite));
        }
        catch(std::exception & e)
        {
            throw e;
        }
    }

    /** \brief Remove elements at a given position. Variants work similar, but erase a range or take iterators
     *
     * \param position int element to erase, starting from 0 to size
     * \return void
     *
     */
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::removeSite(int position)
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
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::removeSite(int start, int end)
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
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::removeSite(VecGenIter start,VecGenIter end)
    {
        try
        {
            m_isSorted=false;
            VecSites.erase(start,end);
        }
        catch(std::exception &e )
        {
            throw e;
        }
    }
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::removeSite(VecGenIter position)
    {
        try
        {
              m_isSorted=false;
            VecSites.erase(position,position);
        }
        catch(std::exception &e)
        {
            throw e;
        }
    }


    /** \brief
     *
     * \param size const int size of the region to generate
     * \param engine std::mt19937& type of engine to use, see manual for details
     * \param sigma const int Standard deviation if we wish to introduct randomness
     * \param ID const std::string Name of the string
     * \return _BASE_ the Element returned
     *
     */
    template <class _BASE_>
    _BASE_ uGenericNGSChrom<_BASE_>::generateRandomSite(const int size, std::mt19937& engine, const int sigma, const std::string ID) const
    {
        try
        {
            uGenericNGSChrom<_BASE_> emptyExcl;
            return generateRandomSite(size,engine,emptyExcl,sigma,ID);
        }
        catch(std::exception &e)
        {
            throw e;
        }
    }

    template <class _BASE_>
    template <class T2>
    _BASE_ uGenericNGSChrom<_BASE_>::generateRandomSite
    (const int size_,std::mt19937& engine,const uGenericNGSChrom<T2> &exclList, const int sigma, const std::string ID) const
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
                        //TODO Is this valid for our basic tags?
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

    template <class _BASE_>
    template <class T2>
    void uGenericNGSChrom<_BASE_>::addNRandomSite
    (const int size, const int n, std::mt19937& engine, const uGenericNGSChrom<T2>& exclList, const int sigma, const std::string ID)
    {
        //Create each tag and add it.
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

    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::addNRandomSite
    (const int size, const int n, std::mt19937& engine, const int sigma, const std::string ID)
    {
        try
        {
            uGenericNGSChrom<_BASE_> emptyExcl;
            addNRandomSite(size,n, engine,emptyExcl,sigma, ID);
        }
        catch(std::exception &e)
        {
            throw e;
        }
    }


    /** \brief Return average element size
     *
     * \return unsigned long long
     *
     */
    template <class _BASE_>
    unsigned long long uGenericNGSChrom<_BASE_>::avgSiteSize() const
    {
        try
        {
            if (this->count() == 0)
                return 0;
            return sumSiteSize()/this->count();
        }
        catch(std::exception & e)
        {
            throw e;
        }
    };

    template <class _BASE_>
    /** \brief return sum of sizes, including overlapping.
     *
     * \return unsigned long long
     *
     */
    unsigned long long uGenericNGSChrom<_BASE_>::sumSiteSize() const
    {
        try
        {
            return accumulateSitesInfo([](unsigned long long partialSum, _BASE_ item) -> unsigned long long
            {
                return partialSum + item.getLenght();
            }, 0ULL);
        }
        catch(std::exception & e)
        {
            throw e;
        }
    }

//return the smallest site size
    template <class _BASE_>
    unsigned long long uGenericNGSChrom<_BASE_>::minSiteSize() const
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
    };

//Largest site size
    template <class _BASE_>
    unsigned long long uGenericNGSChrom<_BASE_>::maxSiteSize() const
    {
        if (this->count() == 0)
            return 0;
        return maxSite(compareLenght)->getLenght();
    };

//Count how many start at same unique positions
    template <class _BASE_>
    long int uGenericNGSChrom<_BASE_>::countUnique() const
    {
        long int count=0;
        int current;
        std::map<int, int >  myUniqueMap;

        typename std::vector<_BASE_>::iterator iterVec;

        for (iterVec = VecSites.begin() ; iterVec!= VecSites.end(); iterVec++)
        {
            current= iterVec.getStart();
            if (myUniqueMap.count(current)==0)
                myUniqueMap.insert( std::pair<int,int>(current,current));
        }
        count=myUniqueMap.size();

        return count;
    };

//Return a vector containing the size of every
    template <class _BASE_>
    std::vector<long long> uGenericNGSChrom<_BASE_>::returnSiteSizes() const
    {
        return computeOnAllSites([] (_BASE_ elem) -> long long {return elem.getLenght();});
    }

//Output functions
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::outputBedFormat(std::ostream& out)
    {

        // applyOnAllSites(bind(bind(&_BASE_::writeBedToOuput), out));
        applyOnAllSites(bind2nd(mem_fun_ref(&_BASE_::writeBedToOuput), out));
    }



    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::printStats(std::ostream& out) const
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
     * \param position const int& Position to evaluate fro
     * \return typename std::vector<_BASE_>::const_iterator
     *
     */
    template <class _BASE_>
    typename std::vector<_BASE_>::const_iterator uGenericNGSChrom<_BASE_>::findPrecedingSite(const int & position) const
    {
        try
        {
            /**< If unsorted, fail */
            if (m_isSorted==false)
                throw ugene_exception_base() <<string_error("findPrecedingSite called on unsorted vector \n") ;
            if ((sortGetStart==nullptr)||(sortGetEnd==nullptr))
                throw ugene_exception_base() <<string_error(" findPrecedingSite called on chrom without appropriate start or end function\n") ;
            auto comp = [&] (const _BASE_ &item1, const int &item2)
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
        catch (ugene_exception_base & e)
        {

#ifdef DEBUG
            std::cerr << "Calling findPrecedingSite on unsorted vector or you did not provide an approriate get function" <<std::endl;
            std::cerr << "sorted status is" << m_isSorted <<std::endl;
            std::cerr << "is Nullprt "<< (sortGetStart==nullptr) <<std::endl;
#endif
            throw e;
        }
    }
//TODO test this and complimentary
    template <class _BASE_>
    typename std::vector<_BASE_>::const_iterator uGenericNGSChrom<_BASE_>::findNextSite(const int & position) const
    {
        try
        {
            /**< If unsorted, fail */
            if ((m_isSorted==false)||(sortGetStart==nullptr)||(sortGetEnd==nullptr))
                throw ugene_exception_base();
            /**< Return true comparitor if item1 smaller then item 2 */
            auto comp = [&] (const int &item1, const _BASE_ &item2)
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
            std::cerr << "Calling findNextSite on unsorted vector or you did not provide an approriate get function" <<std::endl;
            throw;
        }
    }


    /** \brief return the count of a data subset, based on the current sort type. Please see manual for correct usage of this category of function
     *
     * \param start int start interval value
     * \param end int end interval value
     * \param overlap OverlapType OverlapType from Enum
     * \return int Number of elements in our range
     *
     */
    template <class _BASE_>
    int uGenericNGSChrom<_BASE_>::getSubsetCount(int start, int end, OverlapType overlap) const
    {
        try
        {
            auto pos = this->findPrecedingSite(start);
            /**<  If no tag leftwise, we start at beginning*/
            if (pos==this->end())
                pos=this->begin();
            int tagcount=0;
            for (; pos != this->end(); pos++)
            {
                if (sortGetStart(&(*pos))> end)
                    break;
                if (utility::isOverlap(sortGetStart(&(*pos)), sortGetEnd(&(*pos)),start, end,overlap))
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
     * \param start int Start of range
     * \param end int End of range
     * \param overlap OverlapType Type of overlap
     * \return uGenericNGSChrom<_BASE_> Chrom structure containing our element subset
     *
     */
    template <class _BASE_>
    uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::getSubset(int start, int end, OverlapType overlap) const
    {

        try
        {
            uGenericNGSChrom<_BASE_> returnChrom;
            returnChrom.setChr(this->getChr());

            auto pos = this->findPrecedingSite(start);

            /**<  If no tag leftwise, we start at beginning*/
            if (pos==this->end())
                pos=this->begin();

            typename std::vector<_BASE_>::const_iterator iterVec;
            iterVec=VecSites.begin();

            for (; pos != this->end(); pos++)
            {
                if (sortGetStart(&(*pos))> end)
                    break;
                _BASE_ temp;
                if (utility::isOverlap(sortGetStart(&(*pos)), sortGetEnd(&(*pos)),start, end,overlap))
                    returnChrom.addSite(*pos);
            }

            return returnChrom;
        }
        catch(std::exception & e)
        {

            throw e;
        }

    }

    /** \brief return a subSet based on the current comparison value
     *
     * \param start int
     * \param end int
     * \param overlap OverlapType
     * \return uGenericNGSChrom<_BASE_>
     *
     */
    template <class _BASE_>
    uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::removeSubset(int start, int end, OverlapType overlap)
    {

        try
        {
            uGenericNGSChrom<_BASE_> returnChrom;
            returnChrom.setChr(this->getChr());
            std::vector<int> erasePositions;

            auto pos = this->findPrecedingSite(start);

            /**<  If no tag leftwise, we start at beginning*/
            if (pos==this->end())
                pos=this->begin();

            typename std::vector<_BASE_>::iterator iterVec;
            iterVec=VecSites.begin();


            auto delPos=pos;
            for (; pos != this->end(); pos++)
            {
                if (sortGetStart(&(*pos))> end)
                    break;
                /**< When we find a valid element, go back one step and erase th element */
                if (utility::isOverlap(sortGetStart(&(*pos)), sortGetEnd(&(*pos)),start, end,overlap))
                {
                    returnChrom.addSite(*pos);
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
    template <class _BASE_>
    template <class T2>
    uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::getOverlapping(uGenericNGSChrom<T2> &compareChr,OverlapType overlap) const
    {
        try
        {
            uGenericNGSChrom<_BASE_> returnChr;
            for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
            {
                for(auto compit= compareChr.begin(); compit!=compareChr.end(); compit++)
                {
                    if (utility::isOverlap(it->getStart(), it->getEnd(),compit->getStart(),compit->getEnd(),overlap))
                    {
                        returnChr.addSite(*it);
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
    //TODO should this be in here or in chrom?
    template<typename _BASE_>
    uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::getDistinct(std::string chr, int start, int end, OverlapType options)
    {
        uGenericNGSChrom<_BASE_> returnChrom;

        //If you want to use this, you will need to declare a constructur in the parent class of _CHROM_ to manage a _CHROM_<_BASE_> elementa
        //Copy constructor!
        returnChrom= this->getNotOverlapping(start, end);
        return returnChrom;
    }

    /**< Return the elements of A that do not overlap B */
    template <class _BASE_>
    template <class T2>
    uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::getNotOverlapping(uGenericNGSChrom<T2> &compareChr,OverlapType overlap)const
    {
        try
        {
            uGenericNGSChrom<_BASE_> returnChr;
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
                    returnChr.addSite(*it);
                add =true;
            }
            return returnChr;
        }
        catch(std::exception & e)
        {
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
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::divideItemsIntoNBins(int N, SplitType type)
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
     * \param N int : Size of the members to build in BP
     * \param type SplitType : How to manage splits that are not a multiple of N. Possible : STRICT - IGNORE - FILL - ADD
     * \return void
     *
     */
//TODO check this again
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>::divideItemsIntoBinofSize(int N, SplitType type)
    {
        std::vector<_BASE_> newVector;
        std::vector<uGenericNGS> tempVector;
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
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>:: inferChrSize()
    {
        this->maxSite(comparePos)->getEnd();
    }

} // End of namespace NGS
#endif // UFORMATCHROM_H_INCLUDED
