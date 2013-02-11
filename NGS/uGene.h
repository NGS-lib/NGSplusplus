#ifndef UGENE_H_INCLUDED
#define UGENE_H_INCLUDED

/**< A region, is a fairly generic entity */
#include "uFormats.h"
#include "uTags.h"
#include <limits>
#include "utility/utility.h"
namespace NGS
{
class uToken;
class uRegion;
class uParser;
class uBasicNGS;
class uBasicNGSChrom;
class uBasicNGSExperiment;


class uGene : public uGenericNGS<uGene>
{
public:
    uGene();
    uGene(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD);
    uGene(std::string pChr, long long int pStart, long long int pEnd, StrandDir pstrand, float pScore);
    uGene(std::string pChr, long long int pStart, long long int pEnd, float pScore );

    uGene(uTags);
    uGene(uBasicNGS);
    uGene(uRegion);
    uGene(uToken);

    virtual ~uGene();

    enum class featureType: uint_least16_t
    {
        EXON,INTRON, CODING, NCODING, LOOP, PROMOTER, ENHANCER, CUST_1, CUST_2, CUST_3, CUST_4, CUST_5, CUST_6, CUST_7, CUST_8, CUST_9
    };



private:

            class uFeature{

            long long m_start;
            long long m_end;
            featureType m_type;
            std::string m_ID="";
            short int m_offset=0;
        public:

            uFeature(long long pStart, long long pEnd, featureType pType, short int pOffset, std::string pID );
            long long getStart()const{return m_start;};  /**< Return Start of the feature */
            long long getEnd()const{return m_end;};/**< Return End of the feature */
            featureType getType()const{return m_type;};  /**< Return the type of the feature */
            std::string getID()const{return m_ID;}; /**< Return ID of the feature */

            void setType(featureType pType){m_type=pType;}; /**< Set Type of the feature */
            void setID(std::string pID){m_ID = pID;}; /**< Set string ID of the feature */
            void setStart(long long pStart){m_start = pStart;}; /**< Set Start of the feature */
            void setEnd(long long pEnd){m_end = pEnd;};/**< Set End of the feature */
        };
    std::string m_ident=""; /**<  Name of our gene. */
    std::vector<uFeature> m_featureVector;  /**< List of features associated with our gene. Kept as sorted */

public :

    bool isEqual(const uGene & pCompared)const;
    uGene getCopy()const;
    void addFeature(long long, long long, featureType, short int =0, std::string="");

    void removeFeature(std::vector<uFeature>::const_iterator);
    void removeFeature(std::vector<uFeature>::const_iterator,std::vector<uFeature>::const_iterator);
    bool hasFeatureType(featureType);

    int featureCount(){return m_featureVector.size();};
    typename std::vector<uFeature>::const_iterator featureBegin()const;
    typename std::vector<uFeature>::const_iterator featureEnd()const;

};

class uGeneExperiment;
class uGeneChrom :  public uGenericNGSChrom<uGeneChrom,uGene>
{
public:

    uGeneChrom():uGenericNGSChrom() {};
    uGeneChrom(std::string ourChr):uGenericNGSChrom(ourChr)
    {}
    uGeneChrom(std::string ourChr, long long int lenght):uGenericNGSChrom(ourChr,lenght)
    {}
    uGeneChrom(const uGenericNGSChrom<uGeneChrom,uGene>&);
    uGeneChrom& operator=(const uGeneChrom& copFrom);
    uGeneChrom(const uGeneChrom&);
    uGeneChrom(const std::vector<uGene> & copyVec):uGenericNGSChrom(copyVec){};

    uGeneChrom(uBasicNGSChrom);
    uGeneChrom(uTagsChrom);
    uGeneChrom(uRegionChrom);

    /**< End constructor */

    uGeneChrom getCopy()const;

    typename std::vector<uGene>::const_iterator findNextGeneWithFeature(std::string chr, int position)const;
    typename std::vector<uGene>::const_iterator finPrecedingGeneWithFeature(std::string chr, int position)const;

};

class uGeneExperiment: public uGenericNGSExperiment<uGeneExperiment,uGeneChrom, uGene>
{
public:


    uGeneExperiment& operator=(const uGeneExperiment& copFrom)=default;
    uGeneExperiment(const uGeneExperiment&) = default;
    uGeneExperiment()=default;
    uGeneExperiment getCopy() const;

};
} // End of namespace NGS


#endif // UGENE_H_INCLUDED
