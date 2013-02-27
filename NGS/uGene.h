#ifndef UGENE_H_INCLUDED
#define UGENE_H_INCLUDED

/**< A region, is a fairly generic entity */
#include "uFormats.h"
//#include "uTags.h"
#include <limits>
#include "utility/utility.h"
namespace NGS
{
class uToken;
class uRegion;
class uRegionChrom;
class uRegionExperiment;
class uParser;
class uBasicNGS;
class uBasicNGSChrom;
class uBasicNGSExperiment;
class uTags;
class uTagsChrom;
class uTagsExperiment;

//class uGeneChrom;
//class uGeneExperiment;

    enum class featureType: uint_least16_t
    {
        EXON,INTRON, CODING, NCODING, LOOP, PROMOTER, ENHANCER,UTR3,UTR5, MRNA,OPERON,TRNA,
        CUST_1, CUST_2, CUST_3, CUST_4, CUST_5, CUST_6, CUST_7, CUST_8, CUST_9, OTHER
    };

    featureType mapFeature(const std::string &);
    std::string featureStr(const featureType&);

    static const std::map<std::string,featureType> featureMap{ {"EXON",featureType::EXON},{"INTRON",featureType::INTRON},{"CODING",featureType::CODING},{"NCODING",featureType::NCODING  },
    {"LOOP",featureType::LOOP}, {"ENHANCER",featureType::ENHANCER},
    {"UTR3",featureType::UTR3},{"UTR5",featureType::UTR5},{"MRNA",featureType::MRNA},{"OPERON",featureType::OPERON},{"TRNA",featureType::TRNA} };

class uGene : public uGenericNGS<uGene>
{
public:
    uGene(){};
    uGene(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD);
    uGene(std::string pChr, long long int pStart, long long int pEnd, StrandDir pstrand, float pScore);
    uGene(std::string pChr, long long int pStart, long long int pEnd, float pScore );

    uGene(uTags);
    uGene(uBasicNGS);
    uGene(uRegion);
    uGene(uToken);


    ~uGene(){};
private:



    std::string m_class=""; /**<  Class of our type. */
    std::string m_ID=""; /**<  Name for the gene group. */
    std::string m_transcript=""; /**<  If multipled groups have the same name, transcript ID differentiates them. As such, the ID/Transcript Pair must be unique */
    long long int m_BoundaryStart=0; /**< Earliest position of an associated feature */
    long long int m_BoundaryEnd=0; /**< Latest position of an associated feature */

        class uFeature{
            long long m_start; /**< Start position of the feature */
            long long m_end; /**< End position of the feature */
            featureType m_type; /**< Feature Type, strict */
            std::string m_ID=""; /**< ID of the feature */
            std::string m_class="";/**< Class of the feature */
            short int m_offset=0;
            StrandDir m_strand= StrandDir::FORWARD;
        public:

            uFeature(long long pStart, long long pEnd, featureType pType,std::string pID ,std::string pClass, short int pOffset );
            long long getStart()const{return m_start;};  /**< Return Start of the feature */
            long long getEnd()const{return m_end;};/**< Return End of the feature */




            featureType getType()const{return m_type;};  /**< Return the type of the feature */
            std::string getID()const{return m_ID;}; /**< Return ID of the feature */
            std::string getClass()const{return m_class;}; /**< Return ID of the feature */

            void setStrand(StrandDir );
            void setType(featureType pType){m_type=pType;}; /**< Set Type of the feature */
            void setID(std::string pID){m_ID = pID;}; /**< Set string ID of the feature */
            void setClass(std::string pClass){m_class = pClass;}; /**< Set Class ID */
            void setStart(long long pStart){m_start = pStart;}; /**< Set Start of the feature */
            void setEnd(long long pEnd){m_end = pEnd;};/**< Set End of the feature */

            bool operator==(const uFeature &other) const;
            bool operator!=(const uFeature &other) const;

        };
        std::vector<uFeature> m_featureVector;  /**< List of features associated with our gene. Kept as sorted */



public :

    uToken createToken()const;


    std::string getTranscript()const{return m_transcript;}; /**< Return Transcript ID of the gene */
    std::string getID()const{return m_ID;}; /**< Return ID of the gene */
    std::string getClass()const{return m_class;}; /**< Return class of the gene */
    void setID(std::string pID){m_ID = pID;}; /**< Set string ID of the gene */
    void setClass(std::string pClass){m_class = pClass;}; /**< Set Class ID */
    void setTranscript(std::string pTranscript){m_transcript=pTranscript;}; /**< Set Transcript ID */
    long long getBoundaryStart()const{return m_BoundaryStart;};  /**< Return Start of the feature */
    long long getBoundaryEnd()const{return m_BoundaryEnd;};/**< Return End of the feature */





    bool isOverlappingFeature(long long, long long);
    bool isOverlappingFeature(long long, long long, featureType pType);
    bool isEqual(const uGene & pCompared)const;
    uGene getCopy()const;

    void addFeature(long long, long long,StrandDir pStrand, featureType,std::string="", std::string="", short int=0);

    void removeFeature(std::vector<uFeature>::const_iterator);
    void removeFeature(std::vector<uFeature>::const_iterator,std::vector<uFeature>::const_iterator);
    bool hasFeatureType(featureType) const ;

    unsigned int featureCount()const{return m_featureVector.size();};

    typename std::vector<uFeature>::const_iterator featureBegin()const;
    typename std::vector<uFeature>::const_iterator featureEnd()const;


private :

    void inferBoundary();

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

    typename std::vector<uGene>::const_iterator findNextWithFeature(long long pPosition, featureType pType)const;
    typename std::vector<uGene>::const_iterator findPrecedingWithFeature(long long pPosition, featureType pType)const;

    typename std::vector<uGene>::const_iterator findNextFeature(long long pPosition, featureType pType)const;
    typename std::vector<uGene>::const_iterator finPrecedingFeature(long long pPosition, featureType pType)const;

    long long getIDCount(const std::string & pId, const std::string & pTranscript);
    void addData(uToken);
    void addData(const uGene&);

private:


    typename std::vector<uGene>::iterator findGene(const std::string&,const std::string&);

};

class uGeneExperiment: public uGenericNGSExperiment<uGeneExperiment,uGeneChrom, uGene>
{
public:

    void addData(uToken &);
    void addData(const uGeneChrom&);
    uGeneExperiment& operator=(const uGeneExperiment& copFrom)=default;
    uGeneExperiment(const uGeneExperiment&) = default;
    uGeneExperiment()=default;
    uGeneExperiment getCopy() const;


    typename std::vector<uGene>::const_iterator findNextGeneWithFeature(std::string pChr, long long pPosition, featureType pType)const;
    typename std::vector<uGene>::const_iterator findPrecedingGeneWithFeature(std::string pChr,long long pPosition, featureType pType)const;

};
} // End of namespace NGS


#endif // UGENE_H_INCLUDED
