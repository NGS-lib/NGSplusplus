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

/**< When comparing intervals */
enum class OverlapType
{
    OVERLAP_PARTIAL, OVERLAP_COMPLETE, OVERLAP_CENTRE
};
/**< Used when querying a SAM flag */
enum class SamQuery
{
    IS_PAIRED,ALL_ALIGNED_OK, UNMAPPED, NEXT_UNMAPPED, SEQ_REV_STRAND, SEQ_NEXT_REV_STRAND, FIRST_SEG, LAST_SEG, SECOND_ALIGN, FAIL_QUAL, DUPLICATE
};

/**< Used by our parser and others, defined file types we can load/write */
enum class GenomicFileType
{
    BED, SAM
};


namespace utility
{
/**< Validate if A derives from B */

template<typename Child, typename Parent>
//class IsDerivedFrom
//{
//    static_assert(dynamic_cast<Parent*>(static_cast<Child*>(0)) == nullptr, "First class is not derived from the second");
//};
/**< Some conversion functions */
static std::string convertInt(int number);
static int stringToInt(const std::string &ourStr);
static float stringToNumber( const std::string &ourStr );
/**< Template this */
//static inline  std::string clean_WString(const std::string & input_string);
//static inline std::string numberToString(const int number);
//static inline std::string numberToString(const float number);
//static inline std::string numberToString(const double number);


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
static float getSd(const std::vector<float> & ourVec, const float & mean);
static float gaussianSim(float x1, float x2, float sd);
static float getMean(const std::vector<float> & ourVec);
static void debug_string(std::string input_string);


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

//TODO Make more functional. Among other things, offer a "Next Token and return" option;
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

//TODO better error generation
inline static void loadStream(const std::string & filepath,std::ifstream& file)
{
    try
    {
        file.clear();
        file.open(filepath);
        if( !(file))
        {
            throw(10);
        }
    }
    catch(...)
    {
        throw filepath;
    }

}




//TODO Multi-dimentional Gaussian sim
//One dimensional Gaussian
inline static float gaussianSim(float x1, float x2, float sigma)
{

    float returnvalue;
    //Gaussian similarity function
    //TODO abs
    returnvalue= exp(-( pow((x1-x2),2)/(2*pow(sigma,2))));
    return returnvalue;
}

//Utility function courtesy of the internet.
inline static std::string convertInt(int number)
{
    if (number == 0)
        return "0";
    std::string temp="";
    std::string returnvalue="";
    while (number>0)
    {
        temp+=number%10+48;
        number/=10;
    }
    for (unsigned int i=0; i<temp.length(); i++)
        returnvalue+=temp[temp.length()-i-1];
    return returnvalue;
}

//Utility function courtesy of the internet.
inline static std::string numberToString(const int number)
{
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << number;
    return ss.str();
}

inline static std::string numberToString(const float number)
{
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << number;
    return ss.str();
}

inline static std::string numberToString(const double number)
{
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << number;
    return ss.str();
}


//We return a vector with element 0 is Q1, element 1 is median and element 2 is q3
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
};

inline static bool isInInterval(int pos, int start, int end)
{
    if ((pos>=start)&&(pos<=end))
        return true;
    else
        return false;

}

//Is region X1 inside X2
inline static bool isRegionAInsideRegionB( int A1, int A2, int B1, int B2 )
{

    if (isInInterval( A1, B1, B2 )&&(isInInterval(A2, B1, B2)))
        return true;
    else
        return false;
}

//Do these two overlap
inline static bool checkOverlap(int X1, int X2, int Y1, int Y2)
{

    //Start is inside
    if ((X1>=Y1)&&(X1<=Y2))
        return true;

    //Stop is inside
    if ((X2>=Y1)&&(X2<=Y2))
        return true;

    //Y is englobed by X
    if ((X1<=Y1)&&(X2>=Y2))
        return true;

    return false;
}
//Return how many bp overlap
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

        return overlap;
    }
    else
        return overlap;
}
//By default, partial overlap, otherwise as specified.
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


inline static float getMean(const std::vector<float> & ourVec)
{

    float sum = std::accumulate(ourVec.begin(), ourVec.end(), 0.0);

    return (sum/ourVec.size());
}

//Return as float the SD of a list of elements
//Naturally, this is not precise.
inline static float getSd(const std::vector<float> & ourVec, const float & mean)
{
    float sumsq, sd;
    sumsq=0;
    sd=0;
    std::vector<float> deviations;

    for (auto floatit=ourVec.begin(); floatit < ourVec.end(); floatit++ )
    {
        deviations.push_back(*floatit-mean);
    }
    for (auto floatit=deviations.begin(); floatit < deviations.end(); floatit++ )
    {
        *floatit=(pow(*floatit,2));
        sumsq+=*floatit;
    }
    sd = sumsq/(deviations.size()-1);
    sd= sqrt(sd);

    return sd;
}

//Combine a String and int
inline static std::string concatStringInt(std::string ourstring, int ourInt, bool concatstringleft)
{

    std::stringstream sstm;
    if (concatstringleft)
        sstm << ourstring << ourInt;
    else
        sstm <<ourInt <<ourstring;
    return  sstm.str();

}





inline static int stringToInt( const std::string &ourStr )
{

    std::stringstream ss(ourStr);
    int result;
    return ss >> result ? result : 0;
}

inline static float stringToNumber( const std::string &ourStr )
{
    //character array as argument
    std::stringstream ss(ourStr);
    float result;
    return ss >> result ? result : 0;
}

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

static inline void pause_input()
{
    std::cerr << "Pausing, hit anything to continue" <<std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n' );

}



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

static inline void debug_string(std::string input_string)
{
    std::cerr << "String is " <<input_string <<std::endl;

    if (input_string.size())
    {
        char last = input_string.at(input_string.size()-1);


        if (last=='\n')
            std::cerr << "Last character is return " <<input_string <<std::endl;
        if (last=='\t')
            std::cerr << "Last character is tab " <<input_string <<std::endl;
    }
}
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

}



namespace clustering
{

/**<  Average Euclidean distance of all points*/
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
        /**< Patch from Marco to selectively change max */
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


/**<  Courtesy of the internet : http://www.gamedev.net/topic/406316-hausdorff-distance-between-images/ */
//Distance between two set of points, using position in the vector as X
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
    //vecA std::cerr << "maxDistAB = "<< maxDistAB <<std::endl;
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
