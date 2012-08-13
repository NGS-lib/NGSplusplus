#ifndef UFORMATCHROM_H_INCLUDED
#define UFORMATCHROM_H_INCLUDED

#include <climits>
#include <iostream>
#include <map>
#include "utility/utility.h"
#include <algorithm>
#include <parallel/numeric>
#include <functional>

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

public:

    VecGenIter first(){return VecSites.begin();};
    VecGenIter last(){return VecSites.end();};

    void inferChrSize();
    virtual ~uGenericNGSChrom<_BASE_> (){;}


    int count() const
    {
        return VecSites.size();
    };
    int getChromSize() const
    {
        return chromSize;
    };
    void setChromSize(int chromS)
    {
        try
        {
            if (chromS<0)
                throw 10;
            chromSize= chromS;
        }
        catch(int err)
        {
            throw;
        }
    };
    std::string getChr() const
    {
        return chr;
    };
    void setChr(std::string pchr)
    {
        chr = pchr;
    };


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
    uGenericNGSChrom<_BASE_> getSubset(int start, int end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL) const;
    uGenericNGSChrom<_BASE_> removeSubset(int start, int end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
    std::vector<_BASE_> returnVecData()
    {
        return VecSites;
    };

    void outputBedFormat(std::ostream& out);
    void printStats(std::ostream& out) const;

    void divideItemsIntoNBins(int N, SplitType type=SplitType::STRICT);
    void divideItemsIntoBinofSize(int N, SplitType type=SplitType::STRICT);

    template <class T2>
    _BASE_ generateRandomSite(const int size, std::mt19937& engine, const uGenericNGSChrom<T2> &exclList, const int sigma=0, const std::string ID="") const;
    _BASE_ generateRandomSite(const int size, std::mt19937& engine, const int sigma=0, const std::string ID="") const;

    template <class T2>
    void addNRandomSite(const int size, const int n, std::mt19937& engine, const uGenericNGSChrom<T2> &exclList, const int sigma=0, const std::string ID="");
    void addNRandomSite(const int size, const int n, std::mt19937& engine, const int sigma=0, const std::string ID="");
    uGenericNGSChrom<_BASE_> getOverlapping(uGenericNGSChrom &compareExp);
    template <class T2>
    uGenericNGSChrom<_BASE_> getDistinct(uGenericNGSChrom<T2> &compareExp);
    bool addSite(_BASE_ newSite);
    void removeSite(int position);
    int getSubsetCount(int start, int end, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);

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
    template<class UnaryOperation>
    std::vector<_BASE_> transformAndGetVecData(UnaryOperation unary_op)
    {
        std::vector<_BASE_> copyVec(VecSites);
        for_each(begin(copyVec), end(copyVec), unary_op);
        return copyVec;
    }

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
    template<class UnaryOperation>
    auto computeOnAllSites(UnaryOperation unary_op) -> std::vector<decltype(unary_op(_BASE_()))>
    {
        std::vector<decltype(unary_op(_BASE_()))> results;
        results.reserve(VecSites.size());
        transform(begin(VecSites), end(VecSites), std::back_inserter(results), unary_op);
        return results;
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
    template<class UnaryPredicate>
    std::vector<_BASE_> getSpecificSites(UnaryPredicate pred) const
    {
        std::vector<_BASE_> copyVec {};
        copyVec.reserve(VecSites.size());
        copy_if(begin(VecSites), end(VecSites), std::back_inserter(copyVec), pred);
        return copyVec;
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
        return for_each(begin(VecSites), end(VecSites), f);
    }

    template<class UnaryOperation>
    UnaryOperation applyOnAllSites(const UnaryOperation f) const
    {
        return for_each(begin(VecSites), end(VecSites), f);
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
        return __gnu_parallel::accumulate(begin(VecSites), end(VecSites), init, binary_op, __gnu_parallel::sequential_tag());
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
      * \param comp Compare : Binary comparison operation to perform on the sites collection
      * \return void
      */
    template<class Compare>
    void sortSites(Compare comp)
    {

        return std::sort(begin(VecSites), end(VecSites), comp);
    }

    /** \brief Default sort using the start position as a the comparison point
      *
      * \return void
      */
    void sortSites()
    {
        return sortSites(compareStart);
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
        return is_sorted(begin(VecSites), end(VecSites), comp);
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
        return min_element(begin(VecSites), end(VecSites), comp);
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
        return max_element(begin(VecSites), end(VecSites), comp);
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
        return minmax_element(begin(VecSites), end(VecSites), comp);
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
        return count_if(begin(VecSites), end(VecSites), p);
    }
/**< End STL wrappers */
    uGenericNGSChrom() {};
    uGenericNGSChrom(std::string consString):chr(consString){};
    uGenericNGSChrom(std::string consString, long int size);

private :
    static bool compareStart(const _BASE_ &item1, const _BASE_ &item2)
    {
        return item1.getStart() < item2.getStart();
    }

    static bool compareLenght(const _BASE_ &item1, const _BASE_ &item2)
    {
        return item1.getLenght() < item2.getLenght();
    }

    long int chromSize=0;

protected:
    std::vector<_BASE_> VecSites;
    std::string chr;

    std::vector<long long> returnSiteSizes() const;

    unsigned long long avgSiteSize() const;
    unsigned long long minSiteSize() const;
    unsigned long long maxSiteSize() const;
    int countUnique() const;
    unsigned long long sumSiteSize() const;

    int findPrecedingSite(int position, int low, int high);
    // float CorPeaks(uChipSeq otherChip,int extend);
    std::vector<long long> quartilesSize() const;

    template<typename S, typename R> friend class uGenericNGSExperiment;
    // friend class uGenericNGSExperiment;
};

//Constructor with our chrom.
template <class _BASE_>
uGenericNGSChrom<_BASE_>::uGenericNGSChrom(std::string consString, long int size):chr(consString)
{
    try
    {
        setChromSize(size);
    }
    catch(int err)
    {
        throw;
    }
}

template <class _BASE_>
bool uGenericNGSChrom<_BASE_>::addSite(_BASE_ newSite)
{
    //If we are trying to add an entry to the wrong chromosome;

    chr = newSite.getChr();

    VecSites.push_back(newSite);
    return true;
}

//Since we stock in a vector, this is fairly cost heavy
template <class _BASE_>
void uGenericNGSChrom<_BASE_>::removeSite(int position)
{
    if (position>=this->count())
    {
        std::cerr << "Crashing in removeSite() " << position <<" out of range of " <<this->count() << std::endl;
        abort();
    }

    VecSites.erase(VecSites.begin()+(position));
}
template <class _BASE_>
_BASE_ uGenericNGSChrom<_BASE_>::generateRandomSite
(const int size, std::mt19937& engine, const int sigma, const std::string ID) const
{

    uGenericNGSChrom<_BASE_> emptyExcl;
    return generateRandomSite(size,engine,emptyExcl,sigma,ID);
}

//Generator a random tag that can exist in our chrom size, optional seed
//Use time as random generator if no seed specified
//If you want to generate multiple
template <class _BASE_>
template <class T2>
_BASE_ uGenericNGSChrom<_BASE_>::generateRandomSite
(const int size_,std::mt19937& engine,const uGenericNGSChrom<T2> &exclList, const int sigma, const std::string ID) const
{
    //TODO Sanity check here to make sure it is possible to generate the asked for tag.
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

//Create a given number of sites conforming our parameter sizes
template <class _BASE_>
template <class T2>
void uGenericNGSChrom<_BASE_>::addNRandomSite
(const int size, const int n, std::mt19937& engine, const uGenericNGSChrom<T2>& exclList, const int sigma, const std::string ID)
{
    //Create each tag and add it.
    std::string tempID=ID+"-"+getChr()+"-";
    for(int i=0; i<n; i++)
    {
        VecSites.push_back(generateRandomSite(size, engine,exclList,sigma, utility::concatStringInt(tempID, i)));
    }
    chr = VecSites.back().getChr();
}

//Create a given number of sites conforming our parameter sizes
template <class _BASE_>
void uGenericNGSChrom<_BASE_>::addNRandomSite
(const int size, const int n, std::mt19937& engine, const int sigma, const std::string ID)
{

    uGenericNGSChrom<_BASE_> emptyExcl;
    addNRandomSite(size,n, engine,emptyExcl,sigma, ID);
}

//Returns the average size of the sites in our chromosome.
template <class _BASE_>
unsigned long long uGenericNGSChrom<_BASE_>::avgSiteSize() const
{
    if (this->count() == 0)
        return 0;
    return sumSiteSize()/this->count();
};

//Return the sum of the site of all our sides
//This includes overlapping so does not represent coverage
template <class _BASE_>
unsigned long long uGenericNGSChrom<_BASE_>::sumSiteSize() const
{
    return accumulateSitesInfo([](unsigned long long partialSum, _BASE_ item) -> unsigned long long {
                                return partialSum + item.getLenght();
                                }, 0ULL);
}

//return the smallest site size
template <class _BASE_>
unsigned long long uGenericNGSChrom<_BASE_>::minSiteSize() const
{
    if (this->count() == 0)
        return ULONG_MAX;
    return minSite(compareLenght)->getLenght();
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
int uGenericNGSChrom<_BASE_>::countUnique() const
{
    int count=0;
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

//Return a vector containing 3 elements representing 1st, med and 3rd quartil.
template <class _BASE_>
std::vector<long long> uGenericNGSChrom<_BASE_>::quartilesSize() const
{
    std::vector<long long> vectorSizes;
    std::vector<long long> returnVector;
    int q1, med, q3;
    std::vector<int> quartiles;

    typename std::vector<_BASE_>::iterator iterVec;
    vectorSizes = this->returnSiteSizes();
    //quartile positions.
    q1=vectorSizes.size()*0.25;
    med=vectorSizes.size()*0.5;
    q3=vectorSizes.size()*0.75;

    std::nth_element (vectorSizes.begin(), vectorSizes.begin()+q1, vectorSizes.end());
    returnVector.push_back(vectorSizes.at(q1));

    std::nth_element (vectorSizes.begin(), vectorSizes.begin()+med, vectorSizes.end());
    returnVector.push_back(vectorSizes.at(med));

    std::nth_element (vectorSizes.begin(), vectorSizes.begin()+q3, vectorSizes.end());
    returnVector.push_back(vectorSizes.at(q3));

    return returnVector;
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

    quarts = this->quartilesSize();

    out <<"Number of peaks"<< "\t"<< this->count()<<"\n";
    out <<"Average Peak size:"<< "\t"<< this->avgSiteSize()<<"\n";
    out <<"Median size: "<< "\t"<< quarts.at(1)<<"\n";
    out <<"q1 :" << "\t"<< quarts.at(0)<<"\n";
    out <<"q3 :" << "\t"<< quarts.at(2)<<"\n";
    auto minAndMax = minAndMaxSites(compareLenght);
    out <<"Min Peak size:"<< "\t"<< minAndMax.first->getLenght() <<"\n";
    out <<"Max Peak size:"<< "\t"<< minAndMax.second->getLenght() <<"\n";

}

template <class _BASE_>
int uGenericNGSChrom<_BASE_>::getSubsetCount(int start, int end, OverlapType overlap)
{
    int pos=0;
    //Need to mess around with this later, make sure tags at the same position are not being messed with.

    // pos = this->findPrecedingSite(start, 0 , this->count()-1);

    //If no tag leftwise, we start at beginning
    if (pos==-1)
        pos=0;

    typename std::vector<_BASE_>::iterator iterVec;
    iterVec=VecSites.begin();

    int tagcount=0;
    for (iterVec =(iterVec+pos); iterVec != VecSites.end(); ++iterVec)
    {
        if (iterVec->getStart()> end)
            break;


        if (utility::isOverlap(iterVec->getStart(), iterVec->getEnd(),start, end,overlap))
            tagcount++;
    }

    return tagcount;
}

//Find the site ending previous to our region
template <class _BASE_>
int uGenericNGSChrom<_BASE_>::findPrecedingSite(int position, int low, int high)
{
    int  mid = (low + (high - low) / 2);
    if ((high < low)||((mid+1)>high))
        return -1; // not found


    if (VecSites.at(mid).getEnd() > position)
        return findPrecedingSite(position, low, mid-1);
    else
    {

        if ((VecSites.at(mid).getEnd() < position)&&(VecSites.at(mid+1).getEnd()< position))
            return findPrecedingSite(position, mid+1, high);
        else
        {
            return mid;
        }
    }
}


template <class _BASE_>
uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::getSubset(int start, int end, OverlapType overlap) const
{
    uGenericNGSChrom<_BASE_> returnChrom;
    returnChrom.setChr(this->getChr());

    int pos=0;
    //TODO : This need to be validated, something is buggy with our start or end cuttof.
    //  pos = this->findPrecedingSite(start, 0 , this->count()-1);

    //If no tag leftwise, we start at beginning
    if (pos==-1)
        pos=0;

    typename std::vector<_BASE_>::const_iterator iterVec;
    iterVec=VecSites.begin();

    for (iterVec =(iterVec+pos); iterVec != VecSites.end(); ++iterVec)
    {
        if (iterVec->getStart() > end)
            break;
        _BASE_ temp;

        if (utility::isOverlap(iterVec->getStart(), iterVec->getEnd(),start, end,overlap))
            returnChrom.addSite(*iterVec);
    }

    return returnChrom;
}

/**< What does this do? */
template <class _BASE_>
uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::removeSubset(int start, int end, OverlapType overlap)
{
    uGenericNGSChrom<_BASE_> returnChrom;
    returnChrom.setChr(this->getChr());
    std::vector<int> erasePositions;

    int pos=0;
    //Need to mess around with this later, make sure tags at the same position are not being messed with.
    //pos = this->findPrecedingSite(start, 0 , this->count()-1);

    //If no tag leftwise, we start at beginning
    if (pos==-1)
        pos=0;

    typename std::vector<_BASE_>::iterator iterVec;
    iterVec=VecSites.begin();
    for (iterVec =(iterVec+pos); iterVec != VecSites.end(); ++iterVec)
    {

        if (iterVec->getStart()> end)
            break;


        if (utility::isOverlap(iterVec->getStart(), iterVec->getEnd(),start, end,overlap))
        {
            returnChrom.addSite(*iterVec);
            erasePositions.push_back(iterVec - VecSites.begin());
            returnChrom.addSite(*iterVec);
        }
    }

    sort(erasePositions.begin(), erasePositions.end());

    for(int k=(erasePositions.size()-1); k>=0; k--)
    {
        this->removeSite(erasePositions.at(k));
    }

    return returnChrom;
}

//Return the parts of a Chromosome that overlap another given chrom
template <class _BASE_>
uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::getOverlapping(uGenericNGSChrom &compareChr)
{
    uGenericNGSChrom<_BASE_> returnChr;
    for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
    {
        for(auto compit= compareChr.first(); compit!=compareChr.last(); compit++)
        {
            if (utility::isOverlap(it->getStart(), it->getEnd(),compit->getStart(),compit->getEnd()))
            {
                returnChr.addSite(*it);
                break;
            }
        }
    }
    return returnChr;
}

//Return the parts of a Chromosome that overlap another given chrom
template <class _BASE_>
template <class T2>
uGenericNGSChrom<_BASE_> uGenericNGSChrom<_BASE_>::getDistinct(uGenericNGSChrom<T2> &compareExp)
{
    uGenericNGSChrom<_BASE_> returnChr;
    bool add=true;
    std::vector<_BASE_> VecSitesOther = compareExp.returnVecData();
    for(auto it= VecSites.begin(); it!=VecSites.end(); it++)
    {
        for(auto compit= VecSitesOther.begin(); compit!=VecSitesOther.end(); compit++)
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
   try {
    for( _BASE_& x : VecSites){
        tempVector = move(x.divideIntoNBin(N,type));
        newVector.insert( newVector.end(), tempVector.begin(), tempVector.end() );
    }
    VecSites=move(newVector);
   }
    catch(...){
        std::cerr << "Failed in dividerItemsIntoNBins"<<std::endl;
        throw;
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
   try {
    for( _BASE_& x : VecSites){
     //   tempVector =static_cast< vector<_BASE_> >(move(x.divideIntoBinofSize(N,type)));
        tempVector =(move(x.divideIntoBinofSize(N,type)));
        for ( unsigned int i = 0; i < tempVector.size(); ++i ) {
            newVector.push_back(static_cast<_BASE_>( tempVector.at(i) ));
        }

       // newVector.insert( newVector.end(), tempVector.begin(), tempVector.end() );
    }
    VecSites=move(newVector);
   }
    catch(...){
        std::cerr << "Failed in divideItemsIntoBinofSize"<<std::endl;
        throw;
    }
}
    template <class _BASE_>
    void uGenericNGSChrom<_BASE_>:: inferChrSize(){
   // setChrSize( getMax)


    }


#endif // UFORMATCHROM_H_INCLUDED
