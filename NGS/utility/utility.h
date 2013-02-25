#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED

#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <exception>
#include <iostream>
#include <fstream>
#include <stdexcept>

/**< When comparing intervals */
/**< Overlap_Partial is any overlap betwen A and B */
/**< OVERLAP_COMPLETE mean that A is included in B */
/**< OVERLAP_CENTRE means the center of A is overlap by B */
enum class OverlapType
{
    OVERLAP_PARTIAL, OVERLAP_COMPLETE, OVERLAP_CENTRE
};
/**< Used when querying a SAM flag */
enum class SamQuery
{
    IS_PAIRED,ALL_ALIGNED_OK, UNMAPPED, NEXT_UNMAPPED, SEQ_REV_STRAND, SEQ_NEXT_REV_STRAND, FIRST_SEG, LAST_SEG, SECOND_ALIGN, FAIL_QUAL, DUPLICATE
};

namespace utility
{
/**< Debugging functions */
static inline void stringTocerr(const std::string & value);
/**< Overlap functions */
static bool isOverlap(int X1,int X2, int Y1, int Y2, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
static bool isInInterval(const int pos, const int start, const int end);
static bool checkOverlap(const int X1, const int X2, const int Y1, const int Y2);
static bool isRegionAInsideRegionB( int A1, int A2, int B1, int B2 );

/**< Return the quartiles of a given vector */
static std::vector<float> quartilesofVector(std::vector<float> inputVector);
static std::string concatStringInt(std::string ourstring, int ourInt, bool concatstringleft=true);
/**< Minor statistics */
static float getSd( std::vector<float>  ourVec, const float & mean);
static float gaussianSim(float x1, float x2, float sd);
static float getMean(const std::vector<float> & ourVec);
/**< Tokenizer functions */
static void  GetNextToken( std::string& container, size_t& from,const std::string & line );
static void  GetTokens(std::vector<std::string>& tokens, const std::string & line );
/**< Format tests */
static bool is_posnumber(const std::string& s);

/**< Iterator const functions */
template <class Container>
inline static typename Container::iterator to_mutable_iterator(Container& c, typename Container::const_iterator it);


namespace SAM
{
    /** \brief Check and return if a specified sam flag is set from a received sam flag (int)
     *
     * \param flag const int SamFlag
     * \param toQuery const SamQuery Flag we are checking
     * \return bool True if flag is set
     *
     */
    static inline bool querySamFlag(const int flag, const SamQuery toQuery)
    {
        bool query_result;
        switch(toQuery)
        {
        case SamQuery::IS_PAIRED:
            query_result=(flag&0x1);
            break;
        case SamQuery::ALL_ALIGNED_OK:
            query_result=(flag&0x2);
            break;
        case SamQuery::UNMAPPED:
            query_result=(flag&0x4);
            break;
        case SamQuery::NEXT_UNMAPPED:
            query_result=(flag&0x8);
            break;
        case SamQuery::SEQ_REV_STRAND:
            query_result=(flag&0x10);
            break;
        case SamQuery::SEQ_NEXT_REV_STRAND:
            query_result=(flag&0x20);
            break;
        case SamQuery::FIRST_SEG:
            query_result=(flag&0x40);
            break;
        case SamQuery::LAST_SEG:
            query_result=(flag&0x80);
            break;
        case SamQuery::SECOND_ALIGN:
            query_result=(flag&0x100);
            break;
        case SamQuery::FAIL_QUAL:
            query_result=(flag&0x200);
            break;
        case SamQuery::DUPLICATE:
            query_result=(flag&0x400);
            break;
        }
        return query_result;
    }
}


namespace STRING{

        static std::string concatStringListWithSpaces(){return "";};
        /**< Scheme programming in C++. Who woulda thought. */
        template <typename ...Tail>
        static std::string concatStringListWithSpaces(const std::string& curTrack, Tail&&... tail){
        return (curTrack+" "+concatStringListWithSpaces(std::forward<Tail>(tail)...));
    }


}

  template <class Container>
  inline static typename Container::iterator to_mutable_iterator(Container& c, typename Container::const_iterator it)
   {
        return c.begin() + (it - c.begin());
   }

/** \brief Simple tokenizer class that can sometimes be handy to use
 */
class Tokenizer
{
public:
    static const std::string DELIMITERS;
    Tokenizer(const std::string& str);
    Tokenizer(const std::string& str, const std::string& delimiters);
    bool NextToken();
    bool NextToken(const std::string& delimiters);
    const std::string GetToken() const;
    void Reset();
protected:
    size_t m_offset;
    const std::string m_string;
    std::string m_token="";
    std::string m_delimiters;
};

/**< More Tokenizer options */
/**< Courtesy of code from Sourceforce, adapted */
/**< This is mildly more efficient then returning it, as we can reuse the potential memory allocated to the vector */

inline static void  GetTokens(std::vector<std::string>& tokens, const std::string & line )
{
    tokens.clear();
    std:: string buff;

    size_t from = 0;
    while( from < line.length() )
    {
        GetNextToken( buff, from,line );
        if (buff!="")
            tokens.push_back( buff );
    }
}
/**< Hardcoded delimiters are blank space, tab and return (not neewline) */
inline static void  GetNextToken( std::string& container, size_t& from,const std::string & line )
{
    size_t to = from;
    while( from < line.length() && ( line.at(from) == ' ' || line.at(from) == '\t' || line.at(from) == '\r' ) )
        from++;
    to = from + 1;
    while( to < line.length() && line.at(to) != ' ' && line.at(to) != '\t' && line.at(to) != '\r' )
        to++;
// std::cout <<from <<"\t"<<to<<std::endl;
    container = line.substr( from, to - from );
    from = to;
}



/** \brief Wrapper functions, opens a given path with the stream. Checks if valid and returns error if invalid path
 *
 * \param filepath const std::string& Path to the file we want to open
 * \param file std::ifstream& Stream to open file with
 * \return void
 *
 */
inline static void loadStream(const std::string & filepath,std::ifstream& file)
{

    file.clear();
    file.open(filepath);
    if( !(file))
    {
        std::string error("Invalid filepath to stream, failling in LoadStream()\n Path is: "+filepath+"\n");
        throw std::runtime_error(error.c_str());
    }
}
//TODO Multi-dimentional Gaussian sim
//Check
//One dimensional Gaussian
/** \brief Returns the gaussiam sim between two points for a given Standard deviation
 *
 * \param x2 float Point b
 * \param x1 float Point a
 * \param sigma float SD
 * \return float Between 0 and 1
 *
 */
inline static float gaussianSim(float x1, float x2, float sigma)
{
    return exp(-( pow((x1-x2),2)/(2*pow(sigma,2))));
}



/** \brief Checks if a point is inside two other points on a one dimensional axis
 *
 * \param pos int Point to check
 * \param start int Begning of zone
 * \param end int End of zone
 * \return bool True if point is included between start and end
 *
 */
inline static bool isInInterval(int pos, int start, int end)
{
    if ((pos>=start)&&(pos<=end))
        return true;
    else
        return false;
}

/** \brief Is a set of two points  A contained within a different set of two points B on a 1-D axis
 *
 * \param A1 int Point 1 of region A
 * \param A2 int Point 2 of region A
 * \param B1 int Point 1 of region B
 * \param B2 int point 2 of region B
 * \return bool true if A inside B
 *
 */
inline static bool isRegionAInsideRegionB( int A1, int A2, int B1, int B2 )
{
    if (isInInterval( A1, B1, B2 )&&(isInInterval(A2, B1, B2)))
        return true;
    else
        return false;
}


/** \brief Do two points X overlap two points Y
 *
 * \param X1 int Point 1 of X
 * \param X2 int Point 2 of X
 * \param Y1 int Point 1 of Y
 * \param Y2 int Point 2 of Y
 * \return bool true if X overlap Y
 *
 */
inline static bool checkOverlap(int X1, int X2, int Y1, int Y2)
{
    /**<Start is inside  */
    if ((X1>=Y1)&&(X1<=Y2))
        return true;

    /**< Stop is inside */
    if ((X2>=Y1)&&(X2<=Y2))
        return true;

    /**< Y is englobed by X */
    if ((X1<=Y1)&&(X2>=Y2))
        return true;

    return false;
}

/** \brief Count subsection of segemen X that ov erlap segement Y
 *
 * \param X1 int Point 1 of X
 * \param X2 int Point 2 of X
 * \param Y1 int Point 1 of Y
 * \param Y2 int Point 2 of Y
 * \return Size of overlap segment, can be 0
 *
 */
inline static int overlapCount(int X1, int X2, int Y1, int Y2)
{

    int overlap=0;
    int bigX=X1, smallY=Y1;
    if (checkOverlap(X1, X2, Y1, Y2))
    {
        if (X2>X1)
            bigX =X2;

        if (Y2<Y1)
            smallY=Y2;
        overlap =smallY- bigX;
    }

    return overlap;
}
/** \brief Verify if X overlaps Y according to various definitions of "Overlap". Wrapper function
 *
 * \param X1 int Point 1 of X
 * \param X2 int Point 2 of X
 * \param Y1 int Point 1 of Y
 * \param Y2 int Point 2 of Y
 * \param overlap OverlapType What we call an overlap,
 * \return bool True if overlapping
 *
 */
inline static bool isOverlap(int X1, int X2, int Y1, int Y2, OverlapType overlap)
{

    bool result=false;

    switch(overlap)
    {
    case OverlapType::OVERLAP_COMPLETE:
        if (utility::isRegionAInsideRegionB(X1, X2, Y1, Y2))
        {
            result=true;
        }
        break;
    case OverlapType::OVERLAP_PARTIAL:
        if (utility::checkOverlap(X1, X2, Y1, Y2))
        {
            result=true;
        }
        break;

    case OverlapType::OVERLAP_CENTRE:
        int middle = ((X1+ X2)/2);
        if (utility::isInInterval(middle, Y1, Y2))
        {
            result=true;
        }
        break;
    }
    return result;
}


/** \brief Return mean of a vector of floats
 *
 * \param ourVec const std::vector<float>& Vector
 * \return float mean returned
 *
 */
inline static float getMean(const std::vector<float> & ourVec)
{

    float sum = std::accumulate(ourVec.begin(), ourVec.end(), 0.0);
    return (sum/ourVec.size());
}

/** \brief Standard deviations, This is suceptible to overflow/underflow for large values ( or small ones)
 *
 * \param ourVec const std::vector<float>& Vector of values
 * \param mean const float& Mean
 * \return float Returned standard deviation
 *
 */
inline static float getSd(std::vector<float> ourVec, const float & mean)
{
    float sumsq, sd;
    sumsq=0;
    sd=0;

    for (float & floatVal:ourVec )
    {
        floatVal=(floatVal-mean);
    }
    for (auto floatit=ourVec.begin(); floatit < ourVec.end(); floatit++ )
    {
        *floatit=(pow(*floatit,2));
        sumsq+=*floatit;
    }
    sd = sumsq/(ourVec.size()-1);
    sd= sqrt(sd);

    return sd;
}

/** \brief Take a string and int, return a string with the int concatenated left or right
 *
 * \param ourstring std::string original string
 * \param ourInt int int to add
 * \param concatstringleft bool if true, add int to the right, otherwis eleft
 * \return std::string modified string
 *
 */
inline static std::string concatStringInt(std::string ourstring, int ourInt, bool concatstringleft)
{

    std::stringstream sstm;
    if (concatstringleft)
        sstm << ourstring << ourInt;
    else
        sstm <<ourInt <<ourstring;
    return  sstm.str();

}

/** \brief Return a vector containing the first quartile, median and q3 of a vector.
 *
 * \param inputVector std::vector<float> Vector of elements
 * \return std::vector<float> Vector containing 3 values
 *
 */
inline static std::vector<float> quartilesofVector(std::vector<float> inputVector)
{

    std::vector<float> returnVector;
    int q1, med, q3;
    //quartile positions.

    q1=inputVector.size()*0.25;
    med=inputVector.size()*0.5;
    q3=inputVector.size()*0.75;

    std::nth_element (inputVector.begin(), inputVector.begin()+q1, inputVector.end());
    returnVector.push_back(inputVector.at(q1));

    std::nth_element (inputVector.begin(), inputVector.begin()+med, inputVector.end());
    returnVector.push_back(inputVector.at(med));

    std::nth_element (inputVector.begin(), inputVector.begin()+q3, inputVector.end());
    returnVector.push_back(inputVector.at(q3));

    return returnVector;
}


/** \brief Returns quartiles, deprecated , do not use
 *
 */
inline static std::vector<long int> quartilesSize(std::vector<long int> vectorSizes)
{

    std::vector<long int> tempSizes;
    std::vector<long int> returnVector;
    int q1, med, q3;
    std::vector<int> quartiles;

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
}

/** \brief Call to pause program, waits for input
 *
 * \return void
 *
 */
static inline void pause_input()
{
    std::cerr << "Pausing, hit anything to continue" <<std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n' );

}


/** \brief Debug functions, outputs value to standard error.
 *
 * \param value const int& Value to output. Overloaded for float, double and string
 * \return void
 *
 */
static inline void number_to_cerr(const int & value)
{
    std::cerr << value <<std::endl;
}
static inline void number_to_cerr(const float & value)
{
    std::cerr << value <<std::endl;
}
static inline void number_to_cerr(const double & value)
{
    std::cerr << value <<std::endl;
}

static inline void stringTocerr(const std::string & value)
{
    std::cerr << value <<std::endl;
}



/** \brief "Cleans" a window string for use in Unix format. Rarely used
 *
 * \param input_string const std::string& string to clean
 * \return std::string clean string
 *
 */
static inline std::string clean_WString(const std::string & input_string)
{
    size_t cr_idx;
    std::string return_string;
    cr_idx = input_string.find('\r', 0);
    if (cr_idx != std::string::npos)
    {
        return_string=(input_string.substr(0, cr_idx));
    }

    return return_string;
}

static inline bool is_posnumber(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(),
                                      s.end(), [](char c)
    {
        return !std::isdigit(c);
    }) == s.end();
}

static inline int stoi(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    int toReturn;
    if (ss >> toReturn)
    {
        return toReturn;
    }
    throw std::logic_error("utility::stoi -> invalid value in string.");
}

static inline long long stoll(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    long long toReturn;
    if (ss >> toReturn)
    {
        return toReturn;
    }
    throw std::logic_error("utility::stoll -> invalid value in string.");
}

static inline float stof(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    float toReturn;
    if (ss >> toReturn)
    {
        return toReturn;
    }
    throw std::logic_error("utility::stof -> invalid value in string.");
}

static inline std::string to_string(int v)
{
    std::stringstream ss;
    ss << v;
    return ss.str();
}

static inline std::string to_string(float v)
{
    std::stringstream ss;
    ss << v;
    return ss.str();
}

static inline std::string to_string(long v)
{
    std::stringstream ss;
    ss << v;
    return ss.str();
}

static inline std::string to_string(long long v)
{
    std::stringstream ss;
    ss << v;
    return ss.str();
}

static inline std::string to_string(size_t v)
{
    std::stringstream ss;
    ss << v;
    return ss.str();
}
}

/**< Clustering contains odd and ends, distance measures and such */
namespace clustering
{


/** \brief Returns average Euclidean distance of all points of a vector to another. No normalization
 *
 * \param vecA const std::vector<float> Set A
 * \param vecB const std::vector<float> Set B, must be same size as A
 * \return float average Euclidiean distance
 *
 */
static inline float norm_euclidean_dist(const std::vector<float> vecA, const std::vector<float> vecB)
{
    try
    {

        int sizeA = vecA.size(), sizeB = vecB.size();

        if (sizeA!=sizeB)
            throw 10;

        if (sizeA == 0 || sizeB == 0)
        {
            return 0.0;
        }

        float sumErr=0.0;
        for (unsigned int i=0; i<vecA.size(); i++)
        {
            float euclid_Dist = abs(vecA.at(i)-vecB.at(i));
            sumErr+=euclid_Dist;
        }
        sumErr=sumErr/vecA.size();
        return sumErr;
    }
    catch(...)
    {
        std::cerr << "Crashing in sum_squared_error" <<std::endl;
        throw;
    }

}

/***
* dist: Calculate the distance between two set of points.
* Feng, J, X He, ST Mai, and Claudia Plant. 2011.
* â€œA novel similarity measure for fiber clustering using
* longest common subsequence.â€ Proceedings of the 2011: 1-9.
*
* p_ups: Maximum distance for two points to be close.
* p_delta: Window when comparing points between A and B.
*
* return: edit distance between the set of points A and the set of points B.
*
* Thank you to Marc for original code
*/
static inline float align_distance(const std::vector<float> p_A, const std::vector<float> p_B, const float p_ups, const int p_delta, float penalty=0)
{
    int i, j, i1, i2;
    float dDel, dIns, dSub, distance;

    //d is a table with n+1 rows and m+1 columns
    int n=p_A.size()+1;
    int m=p_B.size()+1;

    float* d = new float[n*m];

    for (i=0; i!=n*m; ++i)
    {
        d[i] = 0;
    }
    /**< Find every local maxima at p_delta distance */
    for (i=1; i!=n; ++i)
    {
        /**< Patch  to selectively change max */
        float ups = p_ups; //Donc, si p_ups vaut 0, c'est comme si on l'utilisait pas.

        for (j= std::max(1,i-p_delta); j < std::min(n,i+p_delta); ++j)
            ups = std::max(ups, std::abs(p_A[j-1]));

        for (j= std::max(1,i-p_delta); j < std::min(m,i+p_delta); ++j)
            ups = std::max(ups,std::abs(p_B[j-1]));
        /**< Patchm, reaplce p_ups with ups */
        i1 = i*m;
        i2 = i1-m;
        for (j= std::max(1,i-p_delta); j < std::min(m,i+p_delta); ++j)
        {
            dDel = d[i2+j];
            dIns = d[i1+(j-1)];
            //**< Attempt a local maxima */

            if (ups!=0)
                dSub = d[i2+(j-1)] + ((ups - abs(p_A[i-1] - p_B[j-1]))/ups);
            else
                dSub = d[i2+(j-1)]+1;

            /**< Normal line */
            //dSub = d[i2+(j-1)] + ((p_ups - abs(p_A[i-1] - p_B[j-1]))/p_ups);

            d[i1+j] = std::max(dSub,std::max(dDel,dIns));
        }
    }



    distance = 1.0-(d[std::min(p_A.size(),p_B.size())*(m+1)]/(float)std::max(p_A.size(),p_B.size()));
    //printf(">%f\n", distance);
    if (std::isnan(distance))
    {
        std::ofstream quickstream("errorMat.txt");
        std::cerr << " result is Nan, here is your matrix.." <<std::endl;
        for (i=0; i<n; i++)
        {
            for (j=0; j<m; j++)
            {
                quickstream<<d[i*n+j] << "\t";
            }
            quickstream <<std::endl;
        }

        utility::pause_input();
    }


    delete[] d;


    return distance;
}


/**<  Courtesy of the internet original code: http://www.gamedev.net/topic/406316-hausdorff-distance-between-images/ */

/** \brief Measure the hausdorff distance between two vector of points, points being the Y value and position in the vector the X value
 *
 * \param vecA const std::vector<float>& Set A
 * \param vecB const std::vector<float>& Set B
 * \return float Hausdorf distance
 *
 */
static float hausdorff(const std::vector<float> & vecA, const std::vector<float> & vecB)
{
    int sizeA = vecA.size(), sizeB = vecB.size();

    if (sizeA == 0 || sizeB == 0)
    {
        return 0.0;
    }

    float heightA=*max_element(vecA.begin(),vecA.end());

    float maxDistAB = 0.0;
    for (float i=0.0; i<vecA.size(); i++)
    {
        float minB = 1000000;
        for (float j=0.0; j<vecB.size(); j++)
        {
            /**< the position in our density vector is also our position on the DNA "axis" */
            float tmpDist = sqrt(pow((i-j),2) + pow((vecA.at(i)-vecB.at(j)),2));

            if (tmpDist < minB)
            {
                minB = tmpDist;
            }
        }

        maxDistAB += minB;
    }
    maxDistAB *= 1.0/((float)vecA.size()*heightA);
// Calculate Hausdorff distance d(B,A)
    float maxDistBA = 0;
    for (float i=0; i<vecB.size(); i++)
    {
        float minA = 100000;
        for (float  j=0; j<vecA.size(); j++)
        {
            float tmpDist = sqrt(pow((i-j),2) + pow((vecB.at(i)-vecA.at(j)),2));

            if (tmpDist < minA)
            {
                minA = tmpDist;
            }
        }

        maxDistBA += minA;
    }
    maxDistBA *= 1.0/((float)vecA.size()*heightA);

    auto returnval = std::max(maxDistAB,maxDistBA);
    return returnval;
}

/** \brief Distance between two equaled sized Score vectors
 *
 * \param binA const Vector A containing density
 * \param binB const Vector B containing density
 * \return float List of scores for each bin
 *
 */
static inline float hausdorffTwoRegions(const std::vector<float> &  vectorA, const std::vector<float> & vectorB )
{
    std::vector<float> returnDistance;
    try
    {
        if (vectorA.size()!=vectorA.size())
            throw 10;

        auto result= clustering::hausdorff(vectorA, vectorB);
        return result;
    }
    catch(...)
    {
        std::cerr << "in hauftsmanTwoBins, bins as not same size" <<std::endl;
        throw;
    }

}


}




#endif // UTILITY_H_INCLUDED
