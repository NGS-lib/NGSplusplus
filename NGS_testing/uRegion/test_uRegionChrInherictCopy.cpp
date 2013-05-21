
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

class StandardChromsRegion
{
public:
	StandardChromsRegion()
	{
		/**< m_uRegionChroms[0]: Empty chrom without a name (""). */
       oneChr.setChr("chr1");

       oneChr.addData(uRegion("chr1",100,200));

       manyChr.setChr("chr1");
       manyChr.addData(uRegion("chr1",100,200));
       manyChr.addData(uRegion("chr1",230,300));
       manyChr.addData(uRegion("chr1",120,250));

       emptyChr.setChr("chr1");
	}
	uRegionChrom oneChr;
    uRegionChrom manyChr;
    uRegionChrom emptyChr;
};


TEST(uRegionGENCHR_applyAndGetChrom, NORMAL){

       StandardChromsRegion testChroms;
       auto functOp = [](uRegion & item){   item.extendSite(20);
       };
       uRegionChrom results=testChroms.manyChr.applyAndGetChrom(functOp);
       EXPECT_TRUE(results.getSite(0).isEqual(uRegion("chr1",80,220)));
       EXPECT_TRUE(results.getSite(1).isEqual(uRegion("chr1",210,320)));
       EXPECT_TRUE(results.getSite(2).isEqual(uRegion("chr1",100,270)));
 }
TEST(uRegionGENCHR_applyAndGetChrom, SIDEEFFECT){
       StandardChromsRegion testChroms;
       long int counter=0;
       auto functOp = [&](uRegion & item){   item.extendSite(20);
       counter+=20;
       };
       uRegionChrom results=testChroms.manyChr.applyAndGetChrom(functOp);
       EXPECT_TRUE(results.getSite(0).isEqual(uRegion("chr1",80,220)));
       EXPECT_TRUE(results.getSite(1).isEqual(uRegion("chr1",210,320)));
       EXPECT_TRUE(results.getSite(2).isEqual(uRegion("chr1",100,270)));
       EXPECT_EQ(results.count(),3);
       EXPECT_EQ(counter,60);
 }
TEST(uRegionGENCHR_applyAndGetChrom, EMPTY){
       StandardChromsRegion testChroms;
       auto functOp = [](uRegion & item){   item.extendSite(20);
       };
       uRegionChrom results=testChroms.emptyChr.applyAndGetChrom(functOp);
       EXPECT_EQ(results.count(),0);
 }

/**< Apply on AllVecData */
TEST(uRegionGENCHR_applyAndGetVecData, NORMAL){
        StandardChromsRegion testChroms;
       auto functOp = [&](uRegion & item){   item.extendSite(20);
       };
       auto results=testChroms.manyChr.applyAndGetVecData(functOp);
       EXPECT_TRUE(results.at(0).isEqual(uRegion("chr1",80,220)));
       EXPECT_TRUE(results.at(1).isEqual(uRegion("chr1",210,320)));
       EXPECT_TRUE(results.at(2).isEqual(uRegion("chr1",100,270)));

 }
TEST(uRegionGENCHR_applyAndGetVecData, SIDEEFFECT){
       StandardChromsRegion testChroms;
       long int counter=0;
       auto functOp = [&](uRegion & item){   item.extendSite(20);
       counter+=20;
       };
       auto results=testChroms.manyChr.applyAndGetVecData(functOp);
       EXPECT_TRUE(results.at(0).isEqual(uRegion("chr1",80,220)));
       EXPECT_TRUE(results.at(1).isEqual(uRegion("chr1",210,320)));
       EXPECT_TRUE(results.at(2).isEqual(uRegion("chr1",100,270)));
       EXPECT_EQ((int)results.size(),3);
       EXPECT_EQ(counter,60);
 }
TEST(uRegionGENCHR_applyAndGetVecData, EMPTY){
       StandardChromsRegion testChroms;
       auto functOp = [](uRegion & item){   item.extendSite(20);
       };
       auto results=testChroms.emptyChr.applyAndGetVecData(functOp);
       EXPECT_EQ((int)results.size(),0);
 }

/**<  computeOnAllSites*/
TEST(uRegionGENCHR_computeOnAllSites, NORMAL){
       StandardChromsRegion testChroms;
       auto functOp = [&](uRegion item)->int{  return item.getLength();};
       auto results=testChroms.manyChr.computeOnAllSites(functOp);
       EXPECT_EQ(results.at(0),testChroms.manyChr.getSite(0).getLength());
       EXPECT_EQ(results.at(1),testChroms.manyChr.getSite(1).getLength());
       EXPECT_EQ(results.at(2),testChroms.manyChr.getSite(2).getLength());
       EXPECT_EQ((int)results.size(),testChroms.manyChr.count());

       EXPECT_EQ(std::accumulate(results.begin(), results.end(), 0), (int)testChroms.manyChr.sumSiteSize());
 }
TEST(uRegionGENCHR_computeOnAllSites, EMPTY){
       StandardChromsRegion testChroms;
       auto functOp = [&](uRegion item)->int{  return item.getLength();};
       auto results=testChroms.emptyChr.computeOnAllSites(functOp);
       EXPECT_EQ((int)results.size(),testChroms.emptyChr.count());
 }

/**<  getSpecificSites*/
TEST(uRegionGENCHR_getSpecificSites, NONECOUNTED){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (item.getLength()>2000);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),0);
 }
TEST(uRegionGENCHR_getSpecificSites, SOMECOUNTED){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (item.getLength()>99);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),2);
 }
TEST(uRegionGENCHR_getSpecificSites, ALLCOUNTED){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (item.getLength()>5);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),3);
 }
 TEST(uRegionGENCHR_getSpecificSites, EMPTY){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (true);
       };
       auto results=testChroms.emptyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),0);
 }
/**<  removeSpecificSites
 */
TEST(uRegionGENCHR_removeSpecificSites, NONEREMOVED){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (item.getLength()>2000); };

       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),3);
 }
TEST(uRegionGENCHR_removeSpecificSites, SOMEREMOVED){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (item.getLength()>99);
       };
       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),1);
 }
TEST(uRegionGENCHR_removeSpecificSites, ALLREMOVED){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (item.getLength()>5);
       };
       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),0);
 }
 TEST(uRegionGENCHR_removeSpecificSites, EMPTY){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (false);
       };
       testChroms.emptyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.emptyChr.count(),0);
 }

/**< applyOnAllSites */
TEST(uRegionGENCHR_applyOnAllSites, EMPTY){

       StandardChromsRegion testChroms;
       auto functOp = [](uRegion & item){   item.extendSite(20);
       };
       EXPECT_NO_THROW(testChroms.manyChr.applyOnAllSites(functOp));
 }
TEST(uRegionGENCHR_applyOnAllSites, NORMAL){
        StandardChromsRegion testChroms;
       auto functOp = [](uRegion & item){   item.extendSite(20);
       };
       testChroms.manyChr.applyOnAllSites(functOp);
       EXPECT_TRUE(testChroms.manyChr.getSite(0).isEqual(uRegion("chr1",80,220)));
       EXPECT_TRUE(testChroms.manyChr.getSite(1).isEqual(uRegion("chr1",210,320)));
       EXPECT_TRUE(testChroms.manyChr.getSite(2).isEqual(uRegion("chr1",100,270)));
 }
TEST(uRegionGENCHR_applyOnAllSites, EXCEPTION){

    StandardChromsRegion testChroms;
       auto functOp = [](uRegion & item){   throw param_throw();
       };

      EXPECT_THROW(testChroms.manyChr.applyOnAllSites(functOp),param_throw);

 }
/**< applyOnAllSites */
TEST(uRegionGENCHR_applyOnAllSitesConst, EMPTY){
       const StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){
       };
       EXPECT_NO_THROW(testChroms.manyChr.applyOnAllSites(functOp));
 }
TEST(uRegionGENCHR_applyOnAllSitesConst, NORMAL){
       const StandardChromsRegion testChroms;
       int siteCount=0;
       auto functOp = [&](const uRegion & item){siteCount+=item.getLength();
       };
       testChroms.manyChr.applyOnAllSites(functOp);
       EXPECT_EQ(siteCount,  (int)testChroms.manyChr.sumSiteSize());
 }
TEST(uRegionGENCHR_applyOnAllSitesConst, EXCEPTION){
 StandardChromsRegion testChroms;
       auto functOp = [](uRegion & item){   throw param_throw();
       };
      EXPECT_THROW(testChroms.manyChr.applyOnAllSites(functOp),param_throw);

 }
/**< accumulateSitesInfo */
TEST(uRegionGENCHR_accumulateSitesInfos, EMPTY){
       const StandardChromsRegion testChroms;
       int siteCount=0;
       auto functOp = [&](int siteCounts,const uRegion & item){ return (siteCounts+=item.getLength());
       };
       EXPECT_EQ(testChroms.emptyChr.accumulateSitesInfo(functOp,siteCount),  (int)testChroms.emptyChr.sumSiteSize());
 }
TEST(uRegionGENCHR_accumulateSitesInfo, NORMAL){
       const StandardChromsRegion testChroms;
       int siteCount=0;
       auto functOp = [&](int siteCounts,const uRegion & item){return siteCounts+=item.getLength();
       };
       EXPECT_EQ(testChroms.manyChr.accumulateSitesInfo(functOp,siteCount),  (int)testChroms.manyChr.sumSiteSize());
 }
TEST(uRegionGENCHR_accumulateSitesInfo, EXCEPTION){

       const StandardChromsRegion testChroms;
       auto functOp = [](int siteCounts,const uRegion & item){   throw param_throw();
       return 0;
       };
       int siteCount =0;
       EXPECT_THROW(testChroms.manyChr.accumulateSitesInfo(functOp,siteCount),param_throw);
 }


/**<countSitesWithProperty */
TEST(uRegionGENCHR_countSitesWithProperty, SOMECOUNTED){
      StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item) mutable {  return (item.getLength()>99);
       StandardChromsRegion testChroms;
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),2);
 }

 TEST(uRegionGENCHR_countSitesWithProperty, ALLCOUNTED){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (item.getLength()>2);
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),3);
 }

 TEST(uRegionGENCHR_countSitesWithProperty, NONECOUNTED){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return (item.getLength()>200000);
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),0);
 }
  TEST(uRegionGENCHR_countSitesWithProperty, EMPTY){
       StandardChromsRegion testChroms;
       auto functOp = [](const uRegion & item){  return true;
       };
       EXPECT_EQ(  testChroms.emptyChr.countSitesWithProperty(functOp),0);
 }

/**< minSite(comp ) */
TEST(uRegionGENCHR_minSite, NORMAL)
{
       StandardChromsRegion testChroms;
       auto minItr=  testChroms.manyChr.minSite(uRegionChrom::comparePos);
       EXPECT_EQ(minItr->getStart(),uRegion("chr1",100,200).getStart() );
}
TEST(uRegionGENCHR_minSite, EMPTY)
{
       StandardChromsRegion testChroms;
       auto minItr=  testChroms.emptyChr.minSite(uRegionChrom::comparePos);
       EXPECT_EQ(minItr, testChroms.emptyChr.end() );
}
TEST(uRegionGENCHR_minSite, CUSTOM)
{
       StandardChromsRegion testChroms;
       auto minItr=  testChroms.manyChr.minSite(uRegionChrom::compareLength);
       EXPECT_EQ(minItr->getStart(), uRegion("chr1",230,300).getStart() );
}
/**< maxSite(comp ) */

TEST(uRegionGENCHR_maxSite, NORMAL)
{
       StandardChromsRegion testChroms;
       auto maxItr=  testChroms.manyChr.maxSite(uRegionChrom::comparePos);
       EXPECT_EQ(maxItr->getStart(),uRegion("chr1",230,300).getStart() );
}
TEST(uRegionGENCHR_maxSite, EMPTY)
{

       StandardChromsRegion testChroms;
       auto minItr=  testChroms.emptyChr.maxSite(uRegionChrom::comparePos);
       EXPECT_EQ(minItr, testChroms.emptyChr.end() );
}
TEST(uRegionGENCHR_maxSite, CUSTOM)
{
       StandardChromsRegion testChroms;
       auto minItr=  testChroms.manyChr.maxSite(uRegionChrom::compareLength);
       EXPECT_EQ(minItr->getStart(), uRegion("chr1",120,250).getStart() );
}
/**< MinMaxSites */
TEST(uRegionGENCHR_minAndMaxSites, VALIDCOMPILE){
       StandardChromsRegion testChroms;
        testChroms.manyChr.minAndMaxSites(uRegionChrom::compareLength);
 }

#define CHROMDIVIDESIZE 500
class RegionChromDivide : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {

    uChromTestOverlap.setChr("chr1");
    uChromTestOverlap.addData(uBasicNGS ("chr1", 300, CHROMDIVIDESIZE));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 199));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 299));
  }
uRegionChrom uChromTestOverlap;
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
TEST(uRegionCHR_avgSiteSize, ONESITE){
       StandardChromsRegion ourChroms;
       ourChroms.oneChr.addData(uRegion("chr1",100,200));
       EXPECT_EQ(101,(int)ourChroms.oneChr.avgSiteSize());
 }
TEST(uRegionCHR_avgSiteSize, NOSITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(0,(int)ourChroms.emptyChr.avgSiteSize());
 }
 TEST(uRegionCHR_avgSiteSize, MANYSITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(((101+71+131)/3),(int)ourChroms.manyChr.avgSiteSize());
 }
/**<  */
 TEST(uRegionCHR_minSiteSize, ONESITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(101,(int)ourChroms.oneChr.minSiteSize());
 }
TEST(uRegionCHR_minSiteSize, NOSITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(0,(int)ourChroms.emptyChr.minSiteSize());
 }
 TEST(uRegionCHR_minSiteSize, MANYSITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(71,(int)ourChroms.manyChr.minSiteSize());
 }
 /**<  */
 TEST(uRegionCHR_maxSiteSize, ONESITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(101,(int)ourChroms.oneChr.maxSiteSize() );
 }
TEST(uRegionCHR_maxSiteSizee, NOSITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(0,(int)ourChroms.emptyChr.maxSiteSize());
 }
 TEST(uRegionCHR_maxSiteSize, MANYSITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(131,(int)ourChroms.manyChr.maxSiteSize());
 }
 /**<  */
 TEST(uRegionCHR_sumSiteSize, ONESITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(101,(int)ourChroms.oneChr.sumSiteSize());
 }
TEST(uRegionCHR_sumSiteSize, NOSITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(0,(int)ourChroms.emptyChr.sumSiteSize());
 }
 TEST(uRegionCHR_sumSiteSize, MANYSITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ((101+71+131),(int)ourChroms.manyChr.sumSiteSize());
 }
/**<  */
 TEST(uRegionCHR_inferChrSize, ONESITE){
       StandardChromsRegion ourChroms;
       ourChroms.oneChr.inferChrSize();
       EXPECT_EQ(200,(int)ourChroms.oneChr.getChromSize() );
 }
TEST(uRegionCHR_inferChrSize, NOSITE){
       StandardChromsRegion ourChroms;
       ourChroms.emptyChr.inferChrSize();
       EXPECT_EQ(0,(int)ourChroms.emptyChr.getChromSize() );
 }
 TEST(uRegionCHR_inferChrSize, MANYSITE){
       StandardChromsRegion ourChroms;
       ourChroms.manyChr.inferChrSize();
       EXPECT_EQ(300,(int)ourChroms.manyChr.getChromSize() );
 }
/**<  Count Unique*/
 TEST(uRegionCHR_countUnique, ONESITE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ(1,(int)ourChroms.oneChr.countUnique() );
 }
TEST(uRegionCHR_countUnique, NOSITE){
         StandardChromsRegion ourChroms;
         EXPECT_EQ(0 , (int)ourChroms.emptyChr.countUnique());
 }
 TEST(uRegionCHR_countUnique, MANYSITENOUNIQUE){
       StandardChromsRegion ourChroms;
       EXPECT_EQ( 3, (int)ourChroms.manyChr.countUnique());
 }
TEST(uRegionCHR_countUnique, MANYSITEWITHUNIQUE){
       uRegionChrom newManyChr("chr1");
       newManyChr.addData(uRegion("chr1",100,200));
       newManyChr.addData(uRegion("chr1",100,200));
       newManyChr.addData(uRegion("chr1",150,200));
       newManyChr.addData(uRegion("chr1",100,2500));
       newManyChr.addData(uRegion("chr1",230,300));
       newManyChr.addData(uRegion("chr1",120,250));
       EXPECT_EQ(5 , (int)newManyChr.countUnique());
 }

/**< Testing DivideItemsIntoNBin */
TEST_F(RegionChromDivide, DIVIDEINTONBINFAILCHROM){
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(7));
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(32));
}
TEST_F(RegionChromDivide, DIVIDEINTONBINCHROMIGNORE){
    //std::cout << uChromTestOverlap.count();
    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::IGNORE);
    EXPECT_EQ(uChromTestOverlap.count(),12);
}
TEST_F(RegionChromDivide, DIVIDEINTOBINEXTEND){

    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::EXTEND);
    EXPECT_EQ(uChromTestOverlap.count(),12);
}
TEST_F(RegionChromDivide, DIVIDEINTOBINADD){

    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::ADD);
    EXPECT_EQ(uChromTestOverlap.count(),13);
}
TEST_F(RegionChromDivide, DIVIDECHROMINTONBINOFSIZEFAIL){
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoBinofSize(30));
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoBinofSize(300));
}
TEST_F(RegionChromDivide, DIVIDECHROMINTONBINOFSIZEIGNORE){

    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::IGNORE);
    EXPECT_EQ(uChromTestOverlap.count(),5);
}
TEST_F(RegionChromDivide, DIVIDECHROMINTOBINOFSIZEEXTEND){
    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::EXTEND);
    EXPECT_EQ(uChromTestOverlap.count(),5);
}
TEST_F(RegionChromDivide, DIVIDECHROMINTOBINOFSIZEADD){
    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::ADD);

    EXPECT_EQ(uChromTestOverlap.count(),8);
}

/**<  Random generation test*/
TEST(uRegionCHR_generateRandomSite, NOCHRSIZE){
       StandardChromsRegion ourChroms;
       std::random_device rd;
       std::mt19937 gen(rd());
       EXPECT_THROW(ourChroms.manyChr.generateRandomSite(100,gen,0,"test"),param_throw );
 }
TEST(uRegionCHR_generateRandomSite, LARGETTHENCHRSIZE){
       StandardChromsRegion ourChroms;
       ourChroms.manyChr.setChromSize(500);
       std::random_device rd;
       std::mt19937 gen(rd());
       EXPECT_THROW(ourChroms.manyChr.generateRandomSite(1000,gen,0,"test"),param_throw );
 }

TEST(uRegionCHR_generateRandomSite, MAKEITEM){
       StandardChromsRegion ourChroms;
       ourChroms.manyChr.setChromSize(500);
       std::random_device rd;
       std::mt19937 gen(rd());
       uRegion aItem= ourChroms.manyChr.generateRandomSite(100,gen,0);
       EXPECT_EQ(100,aItem.getLength());
       uRegion aItem2= ourChroms.manyChr.generateRandomSite(150,gen,0);
       EXPECT_EQ(150,aItem2.getLength());
 }
TEST(uRegionCHR_generateRandomSite, MAKEMANYDIFFERENTITEMS){
       StandardChromsRegion ourChroms;
       ourChroms.manyChr.setChromSize(500);

       std::mt19937 gen(432);
       uRegion aItem= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uRegion aItem2= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uRegion aItem3= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uRegion aItem4= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       EXPECT_EQ(aItem2.getLength(),aItem.getLength());
       EXPECT_EQ(aItem3.getLength(),aItem2.getLength());
       EXPECT_EQ(aItem4.getLength(),aItem3.getLength());

       EXPECT_NE(aItem.getStart(),aItem2.getStart());
       EXPECT_NE(aItem2.getStart(),aItem3.getStart());
       EXPECT_NE(aItem3.getStart(),aItem4.getStart());
 }
TEST(uRegionCHR_generateRandomSite, MAKEMANYEXCLUSIONLIST){
      StandardChromsRegion ourChroms;
      uRegionChrom exclusionChrom("chr1");
      exclusionChrom.addData(uRegion("chr1",100,200));
      exclusionChrom.addData(uRegion("chr1",150,500));
      exclusionChrom.addData(uRegion("chr1",800,1000));
      exclusionChrom.setChromSize(1500);
      ourChroms.emptyChr.setChromSize(1500);
      std::mt19937 gen(432);
      for(size_t i =0; i<=2000; i++)
        ourChroms.emptyChr.addData(ourChroms.emptyChr.generateRandomSite(100,gen,exclusionChrom,4,"test"));


      auto chromOverlapping= ourChroms.emptyChr.getOverlapping(exclusionChrom);
      EXPECT_EQ(chromOverlapping.count(),0);
 }

TEST(uRegionCHR_addNRandomSite, MAKECORRECTNUMBER){
       StandardChromsRegion ourChroms;
       ourChroms.emptyChr.setChromSize(1000);
       std::random_device rd;
       std::mt19937 gen(rd());
       ourChroms.emptyChr.addNRandomSite(50,10,gen,0);
       EXPECT_EQ(10,ourChroms.emptyChr.count());
 }
 TEST(uRegionCHR_addNRandomSite, MAKECORRECTNUMBERMANYDIFFERENT){
       StandardChromsRegion ourChroms;
       ourChroms.emptyChr.setChromSize(1000);
       std::random_device rd;
       std::mt19937 gen(4523);
       ourChroms.emptyChr.addNRandomSite(50,10,gen,0);
       EXPECT_EQ(10,ourChroms.emptyChr.count());
 }
 TEST(uRegionCHR_addNRandomSite, MAKEMANYEXCLUSIONLIST){
   StandardChromsRegion ourChroms;
      uRegionChrom exclusionChrom("chr1");
      exclusionChrom.addData(uRegion("chr1",100,200));
      exclusionChrom.addData(uRegion("chr1",150,500));
      exclusionChrom.addData(uRegion("chr1",800,1000));
      exclusionChrom.setChromSize(1500);
      ourChroms.emptyChr.setChromSize(1500);
      std::mt19937 gen(432);
      ourChroms.emptyChr.addNRandomSite(100,2000,gen,exclusionChrom,4);
      EXPECT_EQ(ourChroms.emptyChr.count(),2000);
      auto chromOverlapping= ourChroms.emptyChr.getOverlapping(exclusionChrom);
      EXPECT_EQ(chromOverlapping.count(),0);
 }

/**< Get overlapping */
 TEST(uRegionCHR_getOverlapping, NORMAL){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",50,100));
       uRegionChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(1,retChr.count());
 }
 TEST(uRegionCHR_getOverlapping, FIRSTEMPTY){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",50,100));
       uRegionChrom retChr=ourChroms.emptyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uRegionCHR_getOverlapping, SECONDEMPTY){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       uRegionChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uRegionCHR_getOverlapping, SOMEOVERLAP){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",0,220));
       uRegionChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(2,retChr.count());
 }
  TEST(uRegionCHR_getOverlapping, ALLOVERLAP){
       StandardChromsRegion ourChroms;
       uRegionChrom twoItems("chr1");
       twoItems.addData(uRegion("chr1",0,220));
       twoItems.addData(uRegion("chr1",230,420));
       uRegionChrom retChr=ourChroms.manyChr.getOverlapping(twoItems);
       EXPECT_EQ(3,retChr.count());
 }
 TEST(uRegionCHR_getOverlapping, NOOVERLAP){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",500,1000));
       uRegionChrom retChr=ourChroms.emptyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uRegionCHR_getOverlapping, DIFFERENTCHR){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr2");
       oneItem.addData(uRegion("chr2",0,1000));
       uRegionChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
/**< Get overlapping Count */

 TEST(uRegionCHR_getOverlappingCount, NORMAL){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",50,100));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),1);
 }
 TEST(uRegionCHR_getOverlappingCount, FIRSTEMPTY){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",50,100));
       EXPECT_EQ(ourChroms.emptyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uRegionCHR_getOverlappingCount, SECONDEMPTY){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uRegionCHR_getOverlappingCount, SOMEOVERLAP){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",0,220));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),2);
 }
  TEST(uRegionCHR_getOverlappingCount, ALLOVERLAP){
       StandardChromsRegion ourChroms;
       uRegionChrom twoItems("chr1");
       twoItems.addData(uRegion("chr1",0,220));
       twoItems.addData(uRegion("chr1",230,420));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(twoItems),3);
 }
 TEST(uRegionCHR_getOverlappingCount, NOOVERLAP){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",500,1000));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uRegionCHR_getOverlappingCount, DIFFERENTCHR){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr2");
       oneItem.addData(uRegion("chr2",0,1000));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }

/**< Get not Overlapping */
 TEST(uRegionCHR_getNotOverlapping, NORMAL){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",50,100));
       uRegionChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),2);
 }
 TEST(uRegionCHR_getNotOverlapping, FIRSTEMPTY){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",50,100));
       uRegionChrom retChr=ourChroms.emptyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),0);
 }
  TEST(uRegionCHR_getNotOverlapping, SECONDEMPTY){
       StandardChromsRegion ourChroms;
       uRegionChrom noItems("chr1");
       uRegionChrom retChr=ourChroms.manyChr.getNotOverlapping(noItems);
       EXPECT_EQ(retChr.count(),3);
 }
  TEST(uRegionCHR_getNotOverlapping, SOMEOVERLAP){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",0,220));
       uRegionChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),1);
 }
  TEST(uRegionCHR_getNotOverlapping, NOOVERLAP){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr1");
       oneItem.addData(uRegion("chr1",500,1000));
       uRegionChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),3);
 }
 TEST(uRegionCHR_getNotOverlapping, DIFFERENTCHR){
       StandardChromsRegion ourChroms;
       uRegionChrom oneItem("chr2");
       oneItem.addData(uRegion("chr2",0,1000));
       uRegionChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),3);
 }

/**<  getDistinct*/
TEST(uRegionCHR_getDistinct, EMPTY)
{
     uRegionChrom emptyChr("chr2");
     uRegionChrom testChrom =emptyChr.getDistinct(500, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uRegionCHR_getDistinct, MULTIPLE)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites();
     uRegionChrom testChrom =ourChroms.manyChr.getDistinct(500, 2000);
     EXPECT_EQ(3,testChrom.count());
}
TEST(uRegionCHR_getDistinct, GETNONE)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites();
     uRegionChrom testChrom =ourChroms.manyChr.getDistinct(0, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uRegionCHR_getDistinct, CUSTOMSORT)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
     uRegionChrom testChrom =ourChroms.manyChr.getDistinct(80, 130);
     EXPECT_EQ(2,testChrom.count());
}
TEST(uRegionCHR_getDistinct, FAILNOSORT)
{
     StandardChromsRegion ourChroms;
     EXPECT_THROW(ourChroms.manyChr.getDistinct(80, 130),unsorted_throw );
}
TEST(uRegionCHR_getDistinct, MORECUSTOMSORT)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);

     uRegionChrom testChrom =ourChroms.manyChr.getDistinct(80, 130);
     EXPECT_EQ(2,testChrom.count());
}
/**< Test getSubset, RemoveSubset */
TEST(uRegionCHR_getSubset, EMPTY)
{
     StandardChromsRegion ourChroms;
     uRegionChrom testChrom =ourChroms.emptyChr.getSubset(0, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uRegionCHR_getSubset, MANY)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites();
     uRegionChrom testChrom =ourChroms.manyChr.getSubset(0, 2000);
     EXPECT_EQ(3,testChrom.count());
}
TEST(uRegionCHR_getSubset, UNSORTED)
{
     StandardChromsRegion ourChroms;
     EXPECT_THROW(ourChroms.manyChr.getSubset(80, 130),unsorted_throw );
}
TEST(uRegionCHR_getSubset, CUSTOMEMPTY)
{
     StandardChromsRegion ourChroms;
     ourChroms.emptyChr.sortSites(ourChroms.emptyChr.compareLength,&uRegion::getLength);
     uRegionChrom testChrom =ourChroms.emptyChr.getSubset(80, 130);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uRegionCHR_getSubset, CUSTOMESOME)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
     uRegionChrom testChrom =ourChroms.manyChr.getSubset(80, 130);
     EXPECT_EQ(1,testChrom.count());
}

TEST(uRegionCHR_removeSubset, CUSTOMBORDER)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
     uRegionChrom testChrom =ourChroms.manyChr.removeSubset(80, 130);

     EXPECT_EQ(2,ourChroms.manyChr.count());
     EXPECT_EQ(1,testChrom.count());
}

TEST(uRegionCHR_removeSubset, CUSTOMALL)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
     uRegionChrom testChrom =ourChroms.manyChr.removeSubset(0, 200);
     EXPECT_EQ(3,testChrom.count());
}
TEST(uRegionCHR_removeSubset, EMPTY)
{
     StandardChromsRegion ourChroms;
     uRegionChrom testChrom =ourChroms.emptyChr.removeSubset(0, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uRegionCHR_removeSubset, ALLREMOVED)
{

     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites();
     uRegionChrom testChrom =ourChroms.manyChr.removeSubset(0, 2000);
     EXPECT_EQ(0,ourChroms.manyChr.count());
     EXPECT_EQ(3,testChrom.count());
}
TEST(uRegionCHR_removeSubset, SOMEREMOVED)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites();
     uRegionChrom testChrom =ourChroms.manyChr.removeSubset(0, 150);
     EXPECT_EQ(1,ourChroms.manyChr.count());
     EXPECT_EQ(2,testChrom.count());

}

/**<  getSubsetCount*/

TEST(uRegionCHR_getSubsetCount, NORMAL)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites();
     EXPECT_EQ(3,ourChroms.manyChr.getSubsetCount(0, 2000));
}


TEST(uRegionCHR_getSubsetCount, EMPTY)
{
     StandardChromsRegion ourChroms;
     EXPECT_EQ(0,ourChroms.emptyChr.getSubsetCount(0, 2000));
}

TEST(uRegionCHR_getSubsetCount, MANYCUSTOMSORT)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
     EXPECT_EQ(2,ourChroms.manyChr.getSubsetCount(40, 120));
}
/**< AddData */

TEST(uRegionCHR_addData, VALID)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
     ourChroms.manyChr.addData(uRegion("chr1", 2032,3245));
     EXPECT_FALSE(ourChroms.manyChr.getSortedStatus());
     EXPECT_EQ(4,ourChroms.manyChr.count());

}
TEST(uRegionCHR_addData, INVALIDCHR)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
     EXPECT_THROW(ourChroms.manyChr.addData(uRegion("chr2", 2032,3245)),ugene_exception_base);
}


/**<  GetStartFunct, getEndFunct */
TEST(uRegionCHR_getStartFunct, HASBEENSET)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
     std::function<float(const uRegion*)>  my_funct= &uRegion::getLength;
     EXPECT_NE(nullptr,ourChroms.manyChr.getStartFunct());
}
TEST(uRegionCHR_getStartFunct, NOTBEENSET)
{
     StandardChromsRegion ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getStartFunct());
}

TEST(uRegionCHR_getEndFunct, HASBEENSET)
{
     StandardChromsRegion ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength,&uRegion::getLength);
     std::function<float(const uRegion*)>  my_funct= &uRegion::getLength;
     EXPECT_NE(nullptr,ourChroms.manyChr.getEndFunct());
}

TEST(uRegionCHR_getEndFunct, NOTBEENSET)
{
     StandardChromsRegion ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getEndFunct());
}

TEST(uRegionCHR_getCompFunct, BEENSET)
{
      StandardChromsRegion ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength,&uRegion::getLength);
      EXPECT_NE(nullptr,ourChroms.manyChr.getCompFunct());
}
TEST(uRegionCHR_getCompFunct, NOTBEENSET)
{
   StandardChromsRegion ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getCompFunct());
}

/**<  setChromSize*/

TEST(uRegionCHR_setChromSize, VALID)
{
     StandardChromsRegion ourChroms;
     EXPECT_NO_THROW(ourChroms.emptyChr.setChromSize(20000));
     EXPECT_EQ(20000,ourChroms.emptyChr.getChromSize());
}
TEST(uRegionCHR_setChromSize, INVALID_UNDER0)
{
     StandardChromsRegion ourChroms;
     EXPECT_THROW(ourChroms.emptyChr.setChromSize(-200),param_throw);
}

/**< getSite */
TEST(uRegionCHR_getSite, VALIDREQUEST)
{
     StandardChromsRegion ourChroms;
     uRegion oneItem;
     EXPECT_NO_THROW(oneItem=ourChroms.manyChr.getSite(2));
     EXPECT_TRUE(oneItem.isEqual(uRegion("chr1",120,250)));

}
TEST(uRegionCHR_getSite, INVALID)
{
     StandardChromsRegion ourChroms;
     uRegion oneItem;
     EXPECT_THROW(oneItem=ourChroms.manyChr.getSite(20),param_throw);
}
/**< SortSites */
TEST(uRegionCHR_sortSites, DEFAULT)
{
      StandardChromsRegion ourChroms;
      ourChroms.manyChr.sortSites();
      //Check manually the order
      EXPECT_TRUE(ourChroms.manyChr.getSite(0).isEqual(uRegion("chr1",100,200)));
      EXPECT_TRUE(ourChroms.manyChr.getSite(1).isEqual(uRegion("chr1",120,250)));
      EXPECT_TRUE(ourChroms.manyChr.getSite(2).isEqual(uRegion("chr1",230,300)));
}

TEST(uRegionCHR_sortSites, CUSTOM)
{
      StandardChromsRegion ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength);
      //Check manually the order
      EXPECT_TRUE(ourChroms.manyChr.getSite(1).isEqual(uRegion("chr1",100,200)));
      EXPECT_TRUE(ourChroms.manyChr.getSite(2).isEqual(uRegion("chr1",120,250)));
      EXPECT_TRUE(ourChroms.manyChr.getSite(0).isEqual(uRegion("chr1",230,300)));
}


TEST(uRegionCHR_findNext, STANDARD)
{
      StandardChromsRegion ourChroms;
      ourChroms.manyChr.sortSites();
      auto first=ourChroms.manyChr.findNextSite(195);
      EXPECT_TRUE(first->isEqual(uRegion("chr1",230,300)));
}
TEST(uRegionCHR_findNext, CUSTOM)
{
      StandardChromsRegion ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength,&uRegion::getLength);
      auto first=ourChroms.manyChr.findNextSite(125);
      EXPECT_TRUE(first->isEqual(uRegion("chr1",120,250)));
      first=ourChroms.manyChr.findNextSite(0);
      EXPECT_TRUE(first->isEqual(uRegion("chr1",230,300)));
}

TEST(uRegionCHR_findNext, EMPTY)
{
      StandardChromsRegion ourChroms;
      ourChroms.emptyChr.sortSites();
      auto first=ourChroms.emptyChr.findNextSite(120);
      EXPECT_EQ(first,ourChroms.emptyChr.end());
}

TEST(uRegionCHR_findPrec, STANDARD)
{
      StandardChromsRegion ourChroms;
      ourChroms.manyChr.sortSites();
      auto first=ourChroms.manyChr.findPrecedingSite(195);
      EXPECT_TRUE(first->isEqual(uRegion("chr1",120,250)));
}

TEST(uRegionCHR_findPrec, CUSTOM)
{
      StandardChromsRegion ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uRegion::getLength,&uRegion::getLength);
      auto first=ourChroms.manyChr.findPrecedingSite(125);
      EXPECT_TRUE(first->isEqual(uRegion("chr1",100,200)));
      first=ourChroms.manyChr.findPrecedingSite(400);
      EXPECT_TRUE(first->isEqual(uRegion("chr1",120,250)));
}


TEST(uRegionCHR_findPrec, EMPTY)
{
      StandardChromsRegion ourChroms;
      ourChroms.emptyChr.sortSites();
      auto first=ourChroms.emptyChr.findPrecedingSite(120);
      EXPECT_EQ(first,ourChroms.emptyChr.end());
}
/**< isSorted() */
TEST(uRegionCHR_isSorted, NOT)
{
     StandardChromsRegion ourChroms;
     EXPECT_FALSE(ourChroms.manyChr.isSorted());
}
TEST(uRegionCHR_isSorted, EMPTY)
{
       StandardChromsRegion ourChroms;
       EXPECT_TRUE(ourChroms.emptyChr.isSorted());
}
TEST(uRegionCHR_isSorted, NORMAL)
{
       StandardChromsRegion ourChroms;
       ourChroms.manyChr.sortSites();
       EXPECT_TRUE(ourChroms.emptyChr.isSorted());
       EXPECT_FALSE(ourChroms.manyChr.isSorted(uRegionChrom::compareLength));
}
TEST(uRegionCHR_isSorted, CUSTOM)
{
       StandardChromsRegion ourChroms;
       ourChroms.manyChr.sortSites(uRegionChrom::compareLength);
       EXPECT_TRUE(ourChroms.manyChr.isSorted());
       EXPECT_FALSE(ourChroms.manyChr.isSorted(uRegionChrom::compareStart));
}

TEST(uRegionCHR_RemoveSite, SINGLEITERATOR)
{
       StandardChromsRegion ourChroms;
       ourChroms.manyChr.removeSite(ourChroms.manyChr.begin());
       EXPECT_EQ(2,ourChroms.manyChr.count());
       EXPECT_TRUE(ourChroms.manyChr.begin()->isEqual(uRegion("chr1",230,300)));
        StandardChromsRegion newChroms;
        auto itr =newChroms.manyChr.begin();
        itr++;
        newChroms.manyChr.removeSite(itr);
        EXPECT_EQ(2,newChroms.manyChr.count());
        EXPECT_TRUE(newChroms.manyChr.begin()->isEqual(uRegion("chr1",100,200)));
}
TEST(uRegionCHR_RemoveSite, ITR_RANGE)
{
       StandardChromsRegion ourChroms;
       ourChroms.manyChr.removeSite(ourChroms.manyChr.begin(), ourChroms.manyChr.end());
       EXPECT_EQ(0,ourChroms.manyChr.count());
        StandardChromsRegion newChroms;
        auto itr =newChroms.manyChr.begin();
        auto itr2 =newChroms.manyChr.begin();
        itr2+=2;
        newChroms.manyChr.removeSite(itr,itr2);
        EXPECT_EQ(1,newChroms.manyChr.count());
        EXPECT_TRUE(newChroms.manyChr.begin()->isEqual(uRegion("chr1",120,250)));
}


