// ***************************************************************************
// uGene.h.h (c) 2013
// Alexei Nordell-Markovits : Sherbrooke University
// Charles Joly Beauparlant : Laval University
//
//       This file is part of the NGS++ library.
//
//    The NGS++ library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program (lgpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// ***************************************************************************


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

//GetOverlapFeature

    enum class featureType
    {
        EXON,INTRON, CODING, NCODING, LOOP,LOOP_START,LOOP_END, PROMOTER, ENHANCER,UTR3,UTR5, MRNA,OPERON,TRNA,INTER,INTERCNS,INTRONCNS,
        CUST_1, CUST_2, CUST_3, CUST_4, CUST_5, CUST_6, CUST_7, CUST_8, CUST_9, OTHER, NULLFEATURE
    };

    featureType mapFeature(const std::string &);
    std::string featureString(const featureType&);

    static const std::map<std::string,featureType> featureMap{ {"EXON",featureType::EXON},{"INTRON",featureType::INTRON},{"CODING",featureType::CODING},{"NCODING",featureType::NCODING  },
    {"LOOP",featureType::LOOP},{"LOOP_START",featureType::LOOP_START},{"LOOP_END",featureType::LOOP_END}, {"ENHANCER",featureType::ENHANCER},{"INTER",featureType::INTER},{"INTERCRNS",featureType::INTERCNS},{"INTRONCNS",featureType::INTRONCNS},
    {"UTR3",featureType::UTR3},{"UTR5",featureType::UTR5},{"MRNA",featureType::MRNA},{"OPERON",featureType::OPERON},{"TRNA",featureType::TRNA},{"PROMOTER",featureType::PROMOTER} };


    class uFeature{
        long int m_start; /**< Start position of the feature */
        long int m_end; /**< End position of the feature */
        featureType m_type; /**< Feature Type, strict */
        std::string m_ID=""; /**< ID of the feature */
      //  std::string m_class="";/**< Class of the feature */
        short int m_offset=0;
        StrandDir m_strand= StrandDir::FORWARD;
    public:

        uFeature(long int pStart, long int pEnd,StrandDir, featureType pType,std::string pID , short int pOffset );
        long int getStart()const{return m_start;};  /**< Return Start of the feature */
        long int getEnd()const{return m_end;};/**< Return End of the feature */
        short int getOffSet()const{return m_offset;}; /**< Return offset of the feature*/

        featureType getType()const{return m_type;};  /**< Return the type of the feature */
        std::string getID()const{return m_ID;}; /**< Return ID of the feature */
     //   std::string getClass()const{return m_class;}; /**< Return ID of the feature */


        void setOffset(short int pOffset){m_offset= pOffset;};
        void setStrand(StrandDir pStrandir){m_strand=pStrandir;};
        void setType(featureType pType){m_type=pType;}; /**< Set Type of the feature */
        void setID(std::string pID){m_ID = pID;}; /**< Set string ID of the feature */
     //   void setClass(std::string pClass){m_class = pClass;}; /**< Set Class ID */
        void setStart(long int pStart){m_start = pStart;}; /**< Set Start of the feature */
        void setEnd(long int pEnd){m_end = pEnd;};/**< Set End of the feature */

        bool operator==(const uFeature &other) const;
        bool operator!=(const uFeature &other) const;

    };



class uGene : public uGenericNGS<uGene>
{
public:
    uGene(){};
    uGene(std::string pChr, long int pStart, long int pEnd, StrandDir pStrand=StrandDir::FORWARD);
    uGene(std::string pChr, long int pStart, long int pEnd, StrandDir pstrand, float pScore);
    uGene(std::string pChr, long int pStart, long int pEnd, float pScore );
    uGene(std::string pChr, long int pStart, long int pEnd, StrandDir pStrand, std::string pID, std::string pTranscript="");
    uGene(std::string pChr, long int pStart, long int pEnd, StrandDir pStrand, float pScore, std::string pID, std::string pTranscript="");
    uGene(std::string pChr, long int pStart, long int pEnd, std::string pID, std::string pTranscript="");

    uGene(uTags);
    uGene(uBasicNGS);
    uGene(uRegion);
    uGene(uToken);

    ~uGene(){};
private:



   // std::string m_class=""; /**<  Class of our type. */
    std::string m_ID=""; /**<  Name for the gene group. */
    std::string m_transcript=""; /**<  If multipled groups have the same name, transcript ID differentiates them. As such, the ID/Transcript Pair must be unique */
    long  int m_BoundaryStart=0; /**< Earliest position of an associated feature */
    long  int m_BoundaryEnd=0; /**< Latest position of an associated feature */
    std::vector<uFeature> m_featureVector;  /**< List of features associated with our gene. Kept as sorted */


public :

    uToken createToken()const;

    std::string getTranscript()const{return m_transcript;}; /**< Return Transcript ID of the gene */
    std::string getID()const{return m_ID;}; /**< Return ID of the gene */

    void setID(std::string pID){m_ID = pID;}; /**< Set string ID of the gene */
    void setTranscript(std::string pTranscript){m_transcript=pTranscript;}; /**< Set Transcript ID */
    long int getBoundaryStart()const{return m_BoundaryStart;};  /**< Return Start of the feature */
    long int getBoundaryEnd()const{return m_BoundaryEnd;};/**< Return End of the feature */

    bool isOverlappingFeature(long int, long int);
    bool isOverlappingFeature(long int, long int, featureType pType);
    bool isEqual(const uGene & pCompared)const;
    uGene getCopy()const;

    void addFeature(long int, long int,StrandDir pStrand, featureType, std::string="", short int=0);

    void removeFeature(std::vector<uFeature>::const_iterator);
    void removeFeature(std::vector<uFeature>::const_iterator,std::vector<uFeature>::const_iterator);
    bool hasFeatureType(featureType) const ;

    unsigned long int featureCount(const featureType & pFeature=featureType::NULLFEATURE)const;

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
    uGeneChrom(std::string ourChr, long int length):uGenericNGSChrom(ourChr,length)
    {}
    uGeneChrom(const uGenericNGSChrom<uGeneChrom,uGene>&);
    uGeneChrom& operator=(const uGeneChrom& copFrom);
    uGeneChrom(const uGeneChrom&);

    uGeneChrom(const std::vector<uGene> & copyVec);
    uGeneChrom(const std::string &, const std::vector<uGene> & copyVec);
    uGeneChrom(uBasicNGSChrom);
    uGeneChrom(uTagsChrom);
    uGeneChrom(uRegionChrom);

    /**< End constructor */

    typename std::vector<uGene>::const_iterator findGene(const std::string&,const std::string="")const;
    uGeneChrom getCopy()const;

    unsigned long int featureCount(const featureType &pFeature=featureType::NULLFEATURE)const;

    typename std::vector<uGene>::const_iterator findNextWithFeature(long int pPosition, featureType pType)const;
    //TODO Complete
 //   typename std::vector<uGene>::const_iterator findPrecedingWithFeature(long int pPosition, featureType pType)const;

//TODO complete this
  //  std::pair<typename std::vector<uGene>::const_iterator,typename std::vector<uFeature>::const_iterator> findNextFeature(long int pPosition, featureType pType)const;
 //   std::pair<typename std::vector<uGene>::const_iterator,typename std::vector<uFeature>::const_iterator> finPrecedingFeature(long int pPosition, featureType pType)const;

    long int getIDCount(const std::string & pId, const std::string & pTranscript="")const;
    void addData(const uToken&);
    void addData(const uGene&);

private:

    typename std::vector<uGene>::iterator _findGeneItr(const std::string&,const std::string="");

};

class uGeneExperiment: public uGenericNGSExperiment<uGeneExperiment,uGeneChrom, uGene>
{
public:

    void addData(const uToken&);
    void addData(const uGeneChrom&);
    void addData(const uGene& item){uGenericNGSExperiment::addData(item);};

    uGeneExperiment& operator=(const uGeneExperiment& copFrom)=default;
    uGeneExperiment(const uGeneExperiment&) = default;
    uGeneExperiment()=default;
    uGeneExperiment getCopy() const;

    unsigned long int featureCount(const featureType &pFeature=featureType::NULLFEATURE)const;

    typename std::vector<uGene>::const_iterator findNextWithFeature(std::string pChr, long int pPosition, featureType pType)const;
  //  typename std::vector<uGene>::const_iterator findPrecedingGeneWithFeature(std::string pChr,long int pPosition, featureType pType)const;

};
} // End of namespace NGS


#endif // UGENE_H_INCLUDED
