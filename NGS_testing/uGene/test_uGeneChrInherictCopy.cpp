

/**< Test Common inherited functions */
#include "NGS++.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include <time.h>
#include <gtest/gtest.h>


using namespace std;
using namespace NGS;

class StandardChromsGene
{
public:
	StandardChromsGene()
	{
		/**< m_uGeneChroms[0]: Empty chrom without a name (""). */
       oneChr.setChr("chr1");

       oneChr.addData(uGene("chr1",100,200));

       manyChr.setChr("chr1");
       manyChr.addData(uGene("chr1",100,200));
       manyChr.addData(uGene("chr1",230,300));
       manyChr.addData(uGene("chr1",120,250));

       emptyChr.setChr("chr1");
	}
	uGeneChrom oneChr;
    uGeneChrom manyChr;
    uGeneChrom emptyChr;
};


TEST(uGeneGENCHR_applyAndGetChrom, NORMAL){

       StandardChromsGene testChroms;
       auto functOp = [](uGene & item){   item.extendSite(20);
       };
       uGeneChrom results=testChroms.manyChr.applyAndGetChrom(functOp);
       EXPECT_TRUE(results.getSite(0).isEqual(uGene("chr1",80,220)));
       EXPECT_TRUE(results.getSite(1).isEqual(uGene("chr1",210,320)));
       EXPECT_TRUE(results.getSite(2).isEqual(uGene("chr1",100,270)));
 }
TEST(uGeneGENCHR_applyAndGetChrom, SIDEEFFECT){
       StandardChromsGene testChroms;
       long int counter=0;
       auto functOp = [&](uGene & item){   item.extendSite(20);
       counter+=20;
       };
       uGeneChrom results=testChroms.manyChr.applyAndGetChrom(functOp);
       EXPECT_TRUE(results.getSite(0).isEqual(uGene("chr1",80,220)));
       EXPECT_TRUE(results.getSite(1).isEqual(uGene("chr1",210,320)));
       EXPECT_TRUE(results.getSite(2).isEqual(uGene("chr1",100,270)));
       EXPECT_EQ(results.count(),3);
       EXPECT_EQ(counter,60);
 }
TEST(uGeneGENCHR_applyAndGetChrom, EMPTY){
       StandardChromsGene testChroms;
       auto functOp = [](uGene & item){   item.extendSite(20);
       };
       uGeneChrom results=testChroms.emptyChr.applyAndGetChrom(functOp);
       EXPECT_EQ(results.count(),0);
 }

/**< Apply on AllVecData */
TEST(uGeneGENCHR_applyAndGetVecData, NORMAL){
        StandardChromsGene testChroms;
       auto functOp = [&](uGene & item){   item.extendSite(20);
       };
       auto results=testChroms.manyChr.applyAndGetVecData(functOp);
       EXPECT_TRUE(results.at(0).isEqual(uGene("chr1",80,220)));
       EXPECT_TRUE(results.at(1).isEqual(uGene("chr1",210,320)));
       EXPECT_TRUE(results.at(2).isEqual(uGene("chr1",100,270)));

 }
TEST(uGeneGENCHR_applyAndGetVecData, SIDEEFFECT){
       StandardChromsGene testChroms;
       long int counter=0;
       auto functOp = [&](uGene & item){ item.extendSite(20);
       counter+=20;
       };
       auto results=testChroms.manyChr.applyAndGetVecData(functOp);
       EXPECT_TRUE(results.at(0).isEqual(uGene("chr1",80,220)));
       EXPECT_TRUE(results.at(1).isEqual(uGene("chr1",210,320)));
       EXPECT_TRUE(results.at(2).isEqual(uGene("chr1",100,270)));
       EXPECT_EQ((int)results.size(),3);
       EXPECT_EQ(counter,60);
 }
TEST(uGeneGENCHR_applyAndGetVecData, EMPTY){
       StandardChromsGene testChroms;
       auto functOp = [](uGene & item){   item.extendSite(20);
       };
       auto results=testChroms.emptyChr.applyAndGetVecData(functOp);
       EXPECT_EQ((int)results.size(),0);
 }

/**<  computeOnAllSites*/
TEST(uGeneGENCHR_computeOnAllSites, NORMAL){
       StandardChromsGene testChroms;
       auto functOp = [&](uGene item)->int{  return item.getLength();};
       auto results=testChroms.manyChr.computeOnAllSites(functOp);
       EXPECT_EQ(results.at(0),testChroms.manyChr.getSite(0).getLength());
       EXPECT_EQ(results.at(1),testChroms.manyChr.getSite(1).getLength());
       EXPECT_EQ(results.at(2),testChroms.manyChr.getSite(2).getLength());

       EXPECT_EQ((int)results.size(),testChroms.manyChr.count());


       EXPECT_EQ(std::accumulate(results.begin(), results.end(), 0), (int)testChroms.manyChr.sumSiteSize());
 }
TEST(uGeneGENCHR_computeOnAllSites, EMPTY){
       StandardChromsGene testChroms;
       auto functOp = [&](uGene item)->int{  return item.getLength();};
       auto results=testChroms.emptyChr.computeOnAllSites(functOp);
       EXPECT_EQ((int)results.size(),testChroms.emptyChr.count());
 }

/**<  getSpecificSites*/
TEST(uGeneGENCHR_getSpecificSites, NONECOUNTED){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (item.getLength()>2000);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),0);
 }
TEST(uGeneGENCHR_getSpecificSites, SOMECOUNTED){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (item.getLength()>99);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),2);
 }
TEST(uGeneGENCHR_getSpecificSites, ALLCOUNTED){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (item.getLength()>5);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),3);
 }
 TEST(uGeneGENCHR_getSpecificSites, EMPTY){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (true);
       };
       auto results=testChroms.emptyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),0);
 }
/**<  removeSpecificSites
 */
TEST(uGeneGENCHR_removeSpecificSites, NONEREMOVED){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (item.getLength()>2000); };

       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),3);
 }
TEST(uGeneGENCHR_removeSpecificSites, SOMEREMOVED){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (item.getLength()>99);
       };
       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),1);
 }
TEST(uGeneGENCHR_removeSpecificSites, ALLREMOVED){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (item.getLength()>5);
       };
       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),0);
 }
 TEST(uGeneGENCHR_removeSpecificSites, EMPTY){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (false);
       };
       testChroms.emptyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.emptyChr.count(),0);
 }

/**< applyOnAllSites */
TEST(uGeneGENCHR_applyOnAllSites, EMPTY){

       StandardChromsGene testChroms;
       auto functOp = [](uGene & item){   item.extendSite(20);
       };
       EXPECT_NO_THROW(testChroms.manyChr.applyOnAllSites(functOp));
 }
TEST(uGeneGENCHR_applyOnAllSites, NORMAL){
        StandardChromsGene testChroms;
       auto functOp = [](uGene & item){   item.extendSite(20);
       };
       testChroms.manyChr.applyOnAllSites(functOp);
       EXPECT_TRUE(testChroms.manyChr.getSite(0).isEqual(uGene("chr1",80,220)));
       EXPECT_TRUE(testChroms.manyChr.getSite(1).isEqual(uGene("chr1",210,320)));
       EXPECT_TRUE(testChroms.manyChr.getSite(2).isEqual(uGene("chr1",100,270)));
 }
TEST(uGeneGENCHR_applyOnAllSites, EXCEPTION){

    StandardChromsGene testChroms;
       auto functOp = [](uGene & item){   throw param_throw();
       };

      EXPECT_THROW(testChroms.manyChr.applyOnAllSites(functOp),param_throw);

 }
/**< applyOnAllSites */
TEST(uGeneGENCHR_applyOnAllSitesConst, EMPTY){
       const StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){
       };
       EXPECT_NO_THROW(testChroms.manyChr.applyOnAllSites(functOp));
 }
TEST(uGeneGENCHR_applyOnAllSitesConst, NORMAL){
       const StandardChromsGene testChroms;
       int siteCount=0;
       auto functOp = [&](const uGene & item){siteCount+=item.getLength();
       };
       testChroms.manyChr.applyOnAllSites(functOp);
       EXPECT_EQ(siteCount,  (int)testChroms.manyChr.sumSiteSize());
 }
TEST(uGeneGENCHR_applyOnAllSitesConst, EXCEPTION){
 StandardChromsGene testChroms;
       auto functOp = [](uGene & item){   throw param_throw();
       };
      EXPECT_THROW(testChroms.manyChr.applyOnAllSites(functOp),param_throw);

 }
/**< accumulateSitesInfo */
TEST(uGeneGENCHR_accumulateSitesInfos, EMPTY){
       const StandardChromsGene testChroms;
       int siteCount=0;
       auto functOp = [&](int siteCounts,const uGene & item){ return (siteCounts+=item.getLength());
       };
       EXPECT_EQ((int)testChroms.emptyChr.accumulateSitesInfo(functOp,siteCount),  (int)testChroms.emptyChr.sumSiteSize());
 }
TEST(uGeneGENCHR_accumulateSitesInfo, NORMAL){
       const StandardChromsGene testChroms;
       int siteCount=0;
       auto functOp = [&](int siteCounts,const uGene & item){return siteCounts+=item.getLength();
       };
       EXPECT_EQ((int)testChroms.manyChr.accumulateSitesInfo(functOp,siteCount),  (int)testChroms.manyChr.sumSiteSize());
 }
TEST(uGeneGENCHR_accumulateSitesInfo, EXCEPTION){

       const StandardChromsGene testChroms;
       auto functOp = [](int siteCounts,const uGene & item){   throw param_throw();
       return 0;
       };
       int siteCount =0;
       EXPECT_THROW(testChroms.manyChr.accumulateSitesInfo(functOp,siteCount),param_throw);
 }


/**<countSitesWithProperty */
TEST(uGeneGENCHR_countSitesWithProperty, SOMECOUNTED){
      StandardChromsGene testChroms;
       auto functOp = [](const uGene & item) mutable {  return (item.getLength()>99);
       StandardChromsGene testChroms;
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),2);
 }

 TEST(uGeneGENCHR_countSitesWithProperty, ALLCOUNTED){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (item.getLength()>2);
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),3);
 }

 TEST(uGeneGENCHR_countSitesWithProperty, NONECOUNTED){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return (item.getLength()>200000);
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),0);
 }
  TEST(uGeneGENCHR_countSitesWithProperty, EMPTY){
       StandardChromsGene testChroms;
       auto functOp = [](const uGene & item){  return true;
       };
       EXPECT_EQ(  testChroms.emptyChr.countSitesWithProperty(functOp),0);
 }

/**< minSite(comp ) */
TEST(uGeneGENCHR_minSite, NORMAL)
{
       StandardChromsGene testChroms;
       auto minItr=  testChroms.manyChr.minSite(uGeneChrom::comparePos);
       EXPECT_EQ(minItr->getStart(),uGene("chr1",100,200).getStart() );
}
TEST(uGeneGENCHR_minSite, EMPTY)
{
       StandardChromsGene testChroms;
       auto minItr=  testChroms.emptyChr.minSite(uGeneChrom::comparePos);
       EXPECT_EQ(minItr, testChroms.emptyChr.end() );
}
TEST(uGeneGENCHR_minSite, CUSTOM)
{
       StandardChromsGene testChroms;
       auto minItr=  testChroms.manyChr.minSite(uGeneChrom::compareLength);
       EXPECT_EQ(minItr->getStart(), uGene("chr1",230,300).getStart() );
}
/**< maxSite(comp ) */

TEST(uGeneGENCHR_maxSite, NORMAL)
{
       StandardChromsGene testChroms;
       auto maxItr=  testChroms.manyChr.maxSite(uGeneChrom::comparePos);
       EXPECT_EQ(maxItr->getStart(),uGene("chr1",230,300).getStart() );
}
TEST(uGeneGENCHR_maxSite, EMPTY)
{

       StandardChromsGene testChroms;
       auto minItr=  testChroms.emptyChr.maxSite(uGeneChrom::comparePos);
       EXPECT_EQ(minItr, testChroms.emptyChr.end() );
}
TEST(uGeneGENCHR_maxSite, CUSTOM)
{
       StandardChromsGene testChroms;
       auto minItr=  testChroms.manyChr.maxSite(uGeneChrom::compareLength);
       EXPECT_EQ(minItr->getStart(), uGene("chr1",120,250).getStart() );
}
/**< MinMaxSites */
TEST(uGeneGENCHR_minAndMaxSites, VALIDCOMPILE){
       StandardChromsGene testChroms;
        testChroms.manyChr.minAndMaxSites(uGeneChrom::compareLength);
 }

#define CHROMDIVIDESIZE 500
class GeneChromDivide : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {

    uChromTestOverlap.setChr("chr1");
    uChromTestOverlap.addData(uBasicNGS ("chr1", 300, CHROMDIVIDESIZE));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 199));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 299));
  }
uGeneChrom uChromTestOverlap;
};





/**< Non-Generic */

/*
 * Setters/Getters testing - start/end (positions)
 *	Valid cases:
 *		SETGETSTART
 *		SETGETEND
 *	Invalid cases:
 *		SETSTARTILLEGAL
 *		SETSENDILLEGAL
 */
TEST(uGeneCHR_avgSiteSize, ONESITE){
       StandardChromsGene ourChroms;
       ourChroms.oneChr.addData(uGene("chr1",100,200));
       EXPECT_EQ(101,(int)ourChroms.oneChr.avgSiteSize());
 }
TEST(uGeneCHR_avgSiteSize, NOSITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(0,(int)ourChroms.emptyChr.avgSiteSize());
 }
 TEST(uGeneCHR_avgSiteSize, MANYSITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(((101+71+131)/3),(int)ourChroms.manyChr.avgSiteSize());
 }
/**<  */
 TEST(uGeneCHR_minSiteSize, ONESITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(101,(int)ourChroms.oneChr.minSiteSize());
 }
TEST(uGeneCHR_minSiteSize, NOSITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(0,(int)ourChroms.emptyChr.minSiteSize());
 }
 TEST(uGeneCHR_minSiteSize, MANYSITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(71,(int)ourChroms.manyChr.minSiteSize());
 }
 /**<  */
 TEST(uGeneCHR_maxSiteSize, ONESITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(101,(int)ourChroms.oneChr.maxSiteSize() );
 }
TEST(uGeneCHR_maxSiteSizee, NOSITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(0,(int)ourChroms.emptyChr.maxSiteSize());
 }
 TEST(uGeneCHR_maxSiteSize, MANYSITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(131,(int)ourChroms.manyChr.maxSiteSize());
 }
 /**<  */
 TEST(uGeneCHR_sumSiteSize, ONESITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(101,(int)ourChroms.oneChr.sumSiteSize());
 }
TEST(uGeneCHR_sumSiteSize, NOSITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(0,(int)ourChroms.emptyChr.sumSiteSize());
 }
 TEST(uGeneCHR_sumSiteSize, MANYSITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ((101+71+131),(int)ourChroms.manyChr.sumSiteSize());
 }
/**<  */
 TEST(uGeneCHR_inferChrSize, ONESITE){
       StandardChromsGene ourChroms;
       ourChroms.oneChr.inferChrSize();
       EXPECT_EQ(200,ourChroms.oneChr.getChromSize() );
 }
TEST(uGeneCHR_inferChrSize, NOSITE){
       StandardChromsGene ourChroms;
       ourChroms.emptyChr.inferChrSize();
       EXPECT_EQ(0,ourChroms.emptyChr.getChromSize() );
 }
 TEST(uGeneCHR_inferChrSize, MANYSITE){
       StandardChromsGene ourChroms;
       ourChroms.manyChr.inferChrSize();
       EXPECT_EQ(300,ourChroms.manyChr.getChromSize() );
 }
/**<  Count Unique*/
 TEST(uGeneCHR_countUnique, ONESITE){
       StandardChromsGene ourChroms;
       EXPECT_EQ(1,ourChroms.oneChr.countUnique() );
 }
TEST(uGeneCHR_countUnique, NOSITE){
         StandardChromsGene ourChroms;
         EXPECT_EQ(0 , ourChroms.emptyChr.countUnique());
 }
 TEST(uGeneCHR_countUnique, MANYSITENOUNIQUE){
       StandardChromsGene ourChroms;
       EXPECT_EQ( 3, ourChroms.manyChr.countUnique());
 }
TEST(uGeneCHR_countUnique, MANYSITEWITHUNIQUE){
       uGeneChrom newManyChr("chr1");
       newManyChr.addData(uGene("chr1",100,200));
       newManyChr.addData(uGene("chr1",100,200));
       newManyChr.addData(uGene("chr1",150,200));
       newManyChr.addData(uGene("chr1",100,2500));
       newManyChr.addData(uGene("chr1",230,300));
       newManyChr.addData(uGene("chr1",120,250));
       EXPECT_EQ(5 , newManyChr.countUnique());
 }

/**< Testing DivideItemsIntoNBin */
TEST_F(GeneChromDivide, DIVIDEINTONBINFAILCHROM){
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(7));
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(32));
}
TEST_F(GeneChromDivide, DIVIDEINTONBINCHROMIGNORE){
    //std::cout << uChromTestOverlap.count();
    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::IGNORE);
    EXPECT_EQ(uChromTestOverlap.count(),12);
}
TEST_F(GeneChromDivide, DIVIDEINTOBINEXTEND){

    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::EXTEND);
    EXPECT_EQ(uChromTestOverlap.count(),12);
}
TEST_F(GeneChromDivide, DIVIDEINTOBINADD){

    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::ADD);
    EXPECT_EQ(uChromTestOverlap.count(),13);
}
TEST_F(GeneChromDivide, DIVIDECHROMINTONBINOFSIZEFAIL){
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoBinofSize(30));
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoBinofSize(300));
}
TEST_F(GeneChromDivide, DIVIDECHROMINTONBINOFSIZEIGNORE){

    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::IGNORE);
    EXPECT_EQ(uChromTestOverlap.count(),5);
}
TEST_F(GeneChromDivide, DIVIDECHROMINTOBINOFSIZEEXTEND){
    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::EXTEND);
    EXPECT_EQ(uChromTestOverlap.count(),5);
}
TEST_F(GeneChromDivide, DIVIDECHROMINTOBINOFSIZEADD){
    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::ADD);

    EXPECT_EQ(uChromTestOverlap.count(),8);
}

/**<  Random generation test*/
TEST(uGeneCHR_generateRandomSite, NOCHRSIZE){
       StandardChromsGene ourChroms;
       std::random_device rd;
       std::mt19937 gen(rd());
       EXPECT_THROW(ourChroms.manyChr.generateRandomSite(100,gen,0,"test"),param_throw );
 }
TEST(uGeneCHR_generateRandomSite, LARGETTHENCHRSIZE){
       StandardChromsGene ourChroms;
       ourChroms.manyChr.setChromSize(500);
       std::random_device rd;
       std::mt19937 gen(rd());
       EXPECT_THROW(ourChroms.manyChr.generateRandomSite(1000,gen,0,"test"),param_throw );
 }

TEST(uGeneCHR_generateRandomSite, MAKEITEM){
       StandardChromsGene ourChroms;
       ourChroms.manyChr.setChromSize(500);
       std::random_device rd;
       std::mt19937 gen(rd());
       uGene aItem= ourChroms.manyChr.generateRandomSite(100,gen,0);
       EXPECT_EQ(100,aItem.getLength());
       uGene aItem2= ourChroms.manyChr.generateRandomSite(150,gen,0);
       EXPECT_EQ(150,aItem2.getLength());
 }
TEST(uGeneCHR_generateRandomSite, MAKEMANYDIFFERENTITEMS){
       StandardChromsGene ourChroms;
       ourChroms.manyChr.setChromSize(500);

       std::mt19937 gen(432);
       uGene aItem= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uGene aItem2= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uGene aItem3= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uGene aItem4= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       EXPECT_EQ(aItem2.getLength(),aItem.getLength());
       EXPECT_EQ(aItem3.getLength(),aItem2.getLength());
       EXPECT_EQ(aItem4.getLength(),aItem3.getLength());

       EXPECT_NE(aItem.getStart(),aItem2.getStart());
       EXPECT_NE(aItem2.getStart(),aItem3.getStart());
       EXPECT_NE(aItem3.getStart(),aItem4.getStart());
 }
TEST(uGeneCHR_generateRandomSite, MAKEMANYEXCLUSIONLIST){
      StandardChromsGene ourChroms;
      uGeneChrom exclusionChrom("chr1");
      exclusionChrom.addData(uGene("chr1",100,200));
      exclusionChrom.addData(uGene("chr1",150,500));
      exclusionChrom.addData(uGene("chr1",800,1000));
      exclusionChrom.setChromSize(1500);
      ourChroms.emptyChr.setChromSize(1500);
      std::mt19937 gen(432);
      for(size_t i =0; i<=2000; i++)
        ourChroms.emptyChr.addData(ourChroms.emptyChr.generateRandomSite(100,gen,exclusionChrom,4,"test"));


      auto chromOverlapping= ourChroms.emptyChr.getOverlapping(exclusionChrom);
      EXPECT_EQ(chromOverlapping.count(),0);
 }

TEST(uGeneCHR_addNRandomSite, MAKECORRECTNUMBER){
       StandardChromsGene ourChroms;
       ourChroms.emptyChr.setChromSize(1000);
       std::random_device rd;
       std::mt19937 gen(rd());
       ourChroms.emptyChr.addNRandomSite(50,10,gen,0);
       EXPECT_EQ(10,ourChroms.emptyChr.count());
 }
 TEST(uGeneCHR_addNRandomSite, MAKECORRECTNUMBERMANYDIFFERENT){
       StandardChromsGene ourChroms;
       ourChroms.emptyChr.setChromSize(1000);
       std::random_device rd;
       std::mt19937 gen(4523);
       ourChroms.emptyChr.addNRandomSite(50,10,gen,0);
       EXPECT_EQ(10,ourChroms.emptyChr.count());
 }
 TEST(uGeneCHR_addNRandomSite, MAKEMANYEXCLUSIONLIST){
   StandardChromsGene ourChroms;
      uGeneChrom exclusionChrom("chr1");
      exclusionChrom.addData(uGene("chr1",100,200));
      exclusionChrom.addData(uGene("chr1",150,500));
      exclusionChrom.addData(uGene("chr1",800,1000));
      exclusionChrom.setChromSize(1500);
      ourChroms.emptyChr.setChromSize(1500);
      std::mt19937 gen(432);
      ourChroms.emptyChr.addNRandomSite(100,2000,gen,exclusionChrom,4);
      EXPECT_EQ(ourChroms.emptyChr.count(),2000);
      auto chromOverlapping= ourChroms.emptyChr.getOverlapping(exclusionChrom);
      EXPECT_EQ(chromOverlapping.count(),0);
 }

/**< Get overlapping */
 TEST(uGeneCHR_getOverlapping, NORMAL){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",50,100));
       uGeneChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(1,retChr.count());
 }
 TEST(uGeneCHR_getOverlapping, FIRSTEMPTY){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",50,100));
       uGeneChrom retChr=ourChroms.emptyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uGeneCHR_getOverlapping, SECONDEMPTY){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       uGeneChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uGeneCHR_getOverlapping, SOMEOVERLAP){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",0,220));
       uGeneChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(2,retChr.count());
 }
  TEST(uGeneCHR_getOverlapping, ALLOVERLAP){
       StandardChromsGene ourChroms;
       uGeneChrom twoItems("chr1");
       twoItems.addData(uGene("chr1",0,220));
       twoItems.addData(uGene("chr1",230,420));
       uGeneChrom retChr=ourChroms.manyChr.getOverlapping(twoItems);
       EXPECT_EQ(3,retChr.count());
 }
 TEST(uGeneCHR_getOverlapping, NOOVERLAP){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",500,1000));
       uGeneChrom retChr=ourChroms.emptyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uGeneCHR_getOverlapping, DIFFERENTCHR){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr2");
       oneItem.addData(uGene("chr2",0,1000));
       uGeneChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
/**< Get overlapping Count */

 TEST(uGeneCHR_getOverlappingCount, NORMAL){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",50,100));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),1);
 }
 TEST(uGeneCHR_getOverlappingCount, FIRSTEMPTY){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",50,100));
       EXPECT_EQ(ourChroms.emptyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uGeneCHR_getOverlappingCount, SECONDEMPTY){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uGeneCHR_getOverlappingCount, SOMEOVERLAP){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",0,220));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),2);
 }
  TEST(uGeneCHR_getOverlappingCount, ALLOVERLAP){
       StandardChromsGene ourChroms;
       uGeneChrom twoItems("chr1");
       twoItems.addData(uGene("chr1",0,220));
       twoItems.addData(uGene("chr1",230,420));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(twoItems),3);
 }
 TEST(uGeneCHR_getOverlappingCount, NOOVERLAP){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",500,1000));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uGeneCHR_getOverlappingCount, DIFFERENTCHR){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr2");
       oneItem.addData(uGene("chr2",0,1000));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }

/**< Get not Overlapping */
 TEST(uGeneCHR_getNotOverlapping, NORMAL){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",50,100));
       uGeneChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),2);
 }
 TEST(uGeneCHR_getNotOverlapping, FIRSTEMPTY){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",50,100));
       uGeneChrom retChr=ourChroms.emptyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),0);
 }
  TEST(uGeneCHR_getNotOverlapping, SECONDEMPTY){
       StandardChromsGene ourChroms;
       uGeneChrom noItems("chr1");
       uGeneChrom retChr=ourChroms.manyChr.getNotOverlapping(noItems);
       EXPECT_EQ(retChr.count(),3);
 }
  TEST(uGeneCHR_getNotOverlapping, SOMEOVERLAP){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",0,220));
       uGeneChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),1);
 }
  TEST(uGeneCHR_getNotOverlapping, NOOVERLAP){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr1");
       oneItem.addData(uGene("chr1",500,1000));
       uGeneChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),3);
 }
 TEST(uGeneCHR_getNotOverlapping, DIFFERENTCHR){
       StandardChromsGene ourChroms;
       uGeneChrom oneItem("chr2");
       oneItem.addData(uGene("chr2",0,1000));
       uGeneChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),3);
 }

/**<  getDistinct*/
TEST(uGeneCHR_getDistinct, EMPTY)
{
     uGeneChrom emptyChr("chr2");
     uGeneChrom testChrom =emptyChr.getDistinct(500, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uGeneCHR_getDistinct, MULTIPLE)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites();
     uGeneChrom testChrom =ourChroms.manyChr.getDistinct(500, 2000);
     EXPECT_EQ(3,testChrom.count());
}
TEST(uGeneCHR_getDistinct, GETNONE)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites();
     uGeneChrom testChrom =ourChroms.manyChr.getDistinct(0, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uGeneCHR_getDistinct, CUSTOMSORT)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
     uGeneChrom testChrom =ourChroms.manyChr.getDistinct(80, 130);
     EXPECT_EQ(2,testChrom.count());
}
TEST(uGeneCHR_getDistinct, FAILNOSORT)
{
     StandardChromsGene ourChroms;
     EXPECT_THROW(ourChroms.manyChr.getDistinct(80, 130),unsorted_throw );
}
TEST(uGeneCHR_getDistinct, MORECUSTOMSORT)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);

     uGeneChrom testChrom =ourChroms.manyChr.getDistinct(80, 130);
     EXPECT_EQ(2,testChrom.count());
}
/**< Test getSubset, RemoveSubset */
TEST(uGeneCHR_getSubset, EMPTY)
{
     StandardChromsGene ourChroms;
     uGeneChrom testChrom =ourChroms.emptyChr.getSubset(0, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uGeneCHR_getSubset, MANY)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites();
     uGeneChrom testChrom =ourChroms.manyChr.getSubset(0, 2000);
     EXPECT_EQ(3,testChrom.count());
}
TEST(uGeneCHR_getSubset, UNSORTED)
{
     StandardChromsGene ourChroms;
     EXPECT_THROW(ourChroms.manyChr.getSubset(80, 130),unsorted_throw );
}
TEST(uGeneCHR_getSubset, CUSTOMEMPTY)
{
     StandardChromsGene ourChroms;
     ourChroms.emptyChr.sortSites(ourChroms.emptyChr.compareLength,&uGene::getLength);
     uGeneChrom testChrom =ourChroms.emptyChr.getSubset(80, 130);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uGeneCHR_getSubset, CUSTOMESOME)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
     uGeneChrom testChrom =ourChroms.manyChr.getSubset(80, 130);
     EXPECT_EQ(1,testChrom.count());
}

TEST(uGeneCHR_removeSubset, CUSTOMBORDER)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
     uGeneChrom testChrom =ourChroms.manyChr.removeSubset(80, 130);

     EXPECT_EQ(2,ourChroms.manyChr.count());
     EXPECT_EQ(1,testChrom.count());
}

TEST(uGeneCHR_removeSubset, CUSTOMALL)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
     uGeneChrom testChrom =ourChroms.manyChr.removeSubset(0, 200);
     EXPECT_EQ(3,testChrom.count());
}
TEST(uGeneCHR_removeSubset, EMPTY)
{
     StandardChromsGene ourChroms;
     uGeneChrom testChrom =ourChroms.emptyChr.removeSubset(0, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uGeneCHR_removeSubset, ALLREMOVED)
{

     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites();
     uGeneChrom testChrom =ourChroms.manyChr.removeSubset(0, 2000);
     EXPECT_EQ(0,ourChroms.manyChr.count());
     EXPECT_EQ(3,testChrom.count());
}
TEST(uGeneCHR_removeSubset, SOMEREMOVED)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites();
     uGeneChrom testChrom =ourChroms.manyChr.removeSubset(0, 150);
     EXPECT_EQ(1,ourChroms.manyChr.count());
     EXPECT_EQ(2,testChrom.count());

}

/**<  getSubsetCount*/

TEST(uGeneCHR_getSubsetCount, NORMAL)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites();
     EXPECT_EQ(3,ourChroms.manyChr.getSubsetCount(0, 2000));
}


TEST(uGeneCHR_getSubsetCount, EMPTY)
{
     StandardChromsGene ourChroms;
     EXPECT_EQ(0,ourChroms.emptyChr.getSubsetCount(0, 2000));
}

TEST(uGeneCHR_getSubsetCount, MANYCUSTOMSORT)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
     EXPECT_EQ(2,ourChroms.manyChr.getSubsetCount(40, 120));
}
/**< AddData */

TEST(uGeneCHR_addData, VALID)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
     ourChroms.manyChr.addData(uGene("chr1", 2032,3245));
     EXPECT_FALSE(ourChroms.manyChr.getSortedStatus());
     EXPECT_EQ(4,ourChroms.manyChr.count());

}
TEST(uGeneCHR_addData, INVALIDCHR)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
     EXPECT_THROW(ourChroms.manyChr.addData(uGene("chr2", 2032,3245)),ugene_exception_base);
}


/**<  GetStartFunct, getEndFunct */
TEST(uGeneCHR_getStartFunct, HASBEENSET)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
     std::function<float(const uGene*)>  my_funct= &uGene::getLength;
     EXPECT_NE(nullptr,ourChroms.manyChr.getStartFunct());
}
TEST(uGeneCHR_getStartFunct, NOTBEENSET)
{
     StandardChromsGene ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getStartFunct());
}

TEST(uGeneCHR_getEndFunct, HASBEENSET)
{
     StandardChromsGene ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength,&uGene::getLength);
     std::function<float(const uGene*)>  my_funct= &uGene::getLength;
     EXPECT_NE(nullptr,ourChroms.manyChr.getEndFunct());
}

TEST(uGeneCHR_getEndFunct, NOTBEENSET)
{
     StandardChromsGene ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getEndFunct());
}

TEST(uGeneCHR_getCompFunct, BEENSET)
{
      StandardChromsGene ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength,&uGene::getLength);
      EXPECT_NE(nullptr,ourChroms.manyChr.getCompFunct());
}
TEST(uGeneCHR_getCompFunct, NOTBEENSET)
{
   StandardChromsGene ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getCompFunct());
}

/**<  setChromSize*/

TEST(uGeneCHR_setChromSize, VALID)
{
     StandardChromsGene ourChroms;
     EXPECT_NO_THROW(ourChroms.emptyChr.setChromSize(20000));
     EXPECT_EQ(20000,ourChroms.emptyChr.getChromSize());
}
TEST(uGeneCHR_setChromSize, INVALID_UNDER0)
{
     StandardChromsGene ourChroms;
     EXPECT_THROW(ourChroms.emptyChr.setChromSize(-200),param_throw);
}

/**< getSite */
TEST(uGeneCHR_getSite, VALIDREQUEST)
{
     StandardChromsGene ourChroms;
     uGene oneItem;
     EXPECT_NO_THROW(oneItem=ourChroms.manyChr.getSite(2));
     EXPECT_TRUE(oneItem.isEqual(uGene("chr1",120,250)));

}
TEST(uGeneCHR_getSite, INVALID)
{
     StandardChromsGene ourChroms;
     uGene oneItem;
     EXPECT_THROW(oneItem=ourChroms.manyChr.getSite(20),param_throw);
}
/**< SortSites */
TEST(uGeneCHR_sortSites, DEFAULT)
{
      StandardChromsGene ourChroms;
      ourChroms.manyChr.sortSites();
      //Check manually the order
      EXPECT_TRUE(ourChroms.manyChr.getSite(0).isEqual(uGene("chr1",100,200)));
      EXPECT_TRUE(ourChroms.manyChr.getSite(1).isEqual(uGene("chr1",120,250)));
      EXPECT_TRUE(ourChroms.manyChr.getSite(2).isEqual(uGene("chr1",230,300)));
}

TEST(uGeneCHR_sortSites, CUSTOM)
{
      StandardChromsGene ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength);
      //Check manually the order
      EXPECT_TRUE(ourChroms.manyChr.getSite(1).isEqual(uGene("chr1",100,200)));
      EXPECT_TRUE(ourChroms.manyChr.getSite(2).isEqual(uGene("chr1",120,250)));
      EXPECT_TRUE(ourChroms.manyChr.getSite(0).isEqual(uGene("chr1",230,300)));
}


TEST(uGeneCHR_findNext, STANDARD)
{
      StandardChromsGene ourChroms;
      ourChroms.manyChr.sortSites();
      auto first=ourChroms.manyChr.findNextSite(195);
      EXPECT_TRUE(first->isEqual(uGene("chr1",230,300)));
}
TEST(uGeneCHR_findNext, CUSTOM)
{
      StandardChromsGene ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength,&uGene::getLength);
      auto first=ourChroms.manyChr.findNextSite(125);
      EXPECT_TRUE(first->isEqual(uGene("chr1",120,250)));
      first=ourChroms.manyChr.findNextSite(0);
      EXPECT_TRUE(first->isEqual(uGene("chr1",230,300)));
}

TEST(uGeneCHR_findNext, EMPTY)
{
      StandardChromsGene ourChroms;
      ourChroms.emptyChr.sortSites();
      auto first=ourChroms.emptyChr.findNextSite(120);
      EXPECT_EQ(first,ourChroms.emptyChr.end());
}

TEST(uGeneCHR_findPrec, STANDARD)
{
      StandardChromsGene ourChroms;
      ourChroms.manyChr.sortSites();
      auto first=ourChroms.manyChr.findPrecedingSite(195);
      EXPECT_TRUE(first->isEqual(uGene("chr1",120,250)));
}

TEST(uGeneCHR_findPrec, CUSTOM)
{
      StandardChromsGene ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uGene::getLength,&uGene::getLength);
      auto first=ourChroms.manyChr.findPrecedingSite(125);
      EXPECT_TRUE(first->isEqual(uGene("chr1",100,200)));
      first=ourChroms.manyChr.findPrecedingSite(400);
      EXPECT_TRUE(first->isEqual(uGene("chr1",120,250)));
}


TEST(uGeneCHR_findPrec, EMPTY)
{
      StandardChromsGene ourChroms;
      ourChroms.emptyChr.sortSites();
      auto first=ourChroms.emptyChr.findPrecedingSite(120);
      EXPECT_EQ(first,ourChroms.emptyChr.end());
}
/**< isSorted() */
TEST(uGeneCHR_isSorted, NOT)
{
     StandardChromsGene ourChroms;
     EXPECT_FALSE(ourChroms.manyChr.isSorted());
}
TEST(uGeneCHR_isSorted, EMPTY)
{
       StandardChromsGene ourChroms;
       EXPECT_TRUE(ourChroms.emptyChr.isSorted());
}
TEST(uGeneCHR_isSorted, NORMAL)
{
       StandardChromsGene ourChroms;
       ourChroms.manyChr.sortSites();
       EXPECT_TRUE(ourChroms.emptyChr.isSorted());
       EXPECT_FALSE(ourChroms.manyChr.isSorted(uGeneChrom::compareLength));
}
TEST(uGeneCHR_isSorted, CUSTOM)
{
       StandardChromsGene ourChroms;
       ourChroms.manyChr.sortSites(uGeneChrom::compareLength);
       EXPECT_TRUE(ourChroms.manyChr.isSorted());
       EXPECT_FALSE(ourChroms.manyChr.isSorted(uGeneChrom::compareStart));
}

TEST(uGeneCHR_RemoveSite, SINGLEITERATOR)
{
       StandardChromsGene ourChroms;
       ourChroms.manyChr.removeSite(ourChroms.manyChr.begin());
       EXPECT_EQ(2,ourChroms.manyChr.count());
       EXPECT_TRUE(ourChroms.manyChr.begin()->isEqual(uGene("chr1",230,300)));
        StandardChromsGene newChroms;
        auto itr =newChroms.manyChr.begin();
        itr++;
        newChroms.manyChr.removeSite(itr);
        EXPECT_EQ(2,newChroms.manyChr.count());
        EXPECT_TRUE(newChroms.manyChr.begin()->isEqual(uGene("chr1",100,200)));
}
TEST(uGeneCHR_RemoveSite, ITR_RANGE)
{
       StandardChromsGene ourChroms;
       ourChroms.manyChr.removeSite(ourChroms.manyChr.begin(), ourChroms.manyChr.end());
       EXPECT_EQ(0,ourChroms.manyChr.count());
        StandardChromsGene newChroms;
        auto itr =newChroms.manyChr.begin();
        auto itr2 =newChroms.manyChr.begin();
        itr2+=2;
        newChroms.manyChr.removeSite(itr,itr2);
        EXPECT_EQ(1,newChroms.manyChr.count());
        EXPECT_TRUE(newChroms.manyChr.begin()->isEqual(uGene("chr1",120,250)));
}

