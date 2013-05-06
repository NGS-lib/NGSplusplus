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

#define CHROMDIVIDESIZE 500
using namespace std;
using namespace NGS;

// TODO: printStat (EXPECT_NO_THROW)
// TODO: getSortedStatus
// TODO: get/setChr

#define ONELENGHT 101
class StandardChroms
{
public:
	StandardChroms()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
       oneChr.setChr("chr1");

       oneChr.addData(uBasicNGS("chr1",100,200));

       manyChr.setChr("chr1");
       manyChr.addData(uBasicNGS("chr1",100,200));
       manyChr.addData(uBasicNGS("chr1",230,300));
       manyChr.addData(uBasicNGS("chr1",120,250));

       emptyChr.setChr("chr1");
	}
	uBasicNGSChrom oneChr;
    uBasicNGSChrom manyChr;
    uBasicNGSChrom emptyChr;
};



class ChromDivide : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {

    uChromTestOverlap.setChr("chr1");
    uChromTestOverlap.addData(uBasicNGS ("chr1", 300, CHROMDIVIDESIZE));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 199));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 299));
  }
uBasicNGSChrom uChromTestOverlap;
};



/*
 * Setters/Getters testing - start/end (positions)
 *	Valid cases:
 *		SETGETSTART
 *		SETGETEND
 *	Invalid cases:
 *		SETSTARTILLEGAL
 *		SETSENDILLEGAL
 */
TEST(uBasicNGSCHR_avgSiteSize, ONESITE){
       StandardChroms ourChroms;
       ourChroms.oneChr.addData(uBasicNGS("chr1",100,200));
       EXPECT_EQ(101,ourChroms.oneChr.avgSiteSize());
 }
TEST(uBasicNGSCHR_avgSiteSize, NOSITE){
       StandardChroms ourChroms;
       EXPECT_EQ(0,ourChroms.emptyChr.avgSiteSize());
 }
 TEST(uBasicNGSCHR_avgSiteSize, MANYSITE){
       StandardChroms ourChroms;
       EXPECT_EQ(((101+71+131)/3),ourChroms.manyChr.avgSiteSize());
 }
/**<  */
 TEST(uBasicNGSCHR_minSiteSize, ONESITE){
       StandardChroms ourChroms;
       EXPECT_EQ(101,ourChroms.oneChr.minSiteSize());
 }
TEST(uBasicNGSCHR_minSiteSize, NOSITE){
       StandardChroms ourChroms;
       EXPECT_EQ(0,ourChroms.emptyChr.minSiteSize());
 }
 TEST(uBasicNGSCHR_minSiteSize, MANYSITE){
       StandardChroms ourChroms;
       EXPECT_EQ(71,ourChroms.manyChr.minSiteSize());
 }
 /**<  */
 TEST(uBasicNGSCHR_maxSiteSize, ONESITE){
       StandardChroms ourChroms;
       EXPECT_EQ(101,ourChroms.oneChr.maxSiteSize() );
 }
TEST(uBasicNGSCHR_maxSiteSizee, NOSITE){
       StandardChroms ourChroms;
       EXPECT_EQ(0,ourChroms.emptyChr.maxSiteSize());
 }
 TEST(uBasicNGSCHR_maxSiteSize, MANYSITE){
       StandardChroms ourChroms;
       EXPECT_EQ(131,ourChroms.manyChr.maxSiteSize());
 }
 /**<  */
 TEST(uBasicNGSCHR_sumSiteSize, ONESITE){
       StandardChroms ourChroms;
       EXPECT_EQ(101,ourChroms.oneChr.sumSiteSize());
 }
TEST(uBasicNGSCHR_sumSiteSize, NOSITE){
       StandardChroms ourChroms;
       EXPECT_EQ(0,ourChroms.emptyChr.sumSiteSize());
 }
 TEST(uBasicNGSCHR_sumSiteSize, MANYSITE){
       StandardChroms ourChroms;
       EXPECT_EQ((101+71+131),ourChroms.manyChr.sumSiteSize());
 }
/**<  */
 TEST(uBasicNGSCHR_inferChrSize, ONESITE){
       StandardChroms ourChroms;
       ourChroms.oneChr.inferChrSize();
       EXPECT_EQ(200,ourChroms.oneChr.getChromSize() );
 }
TEST(uBasicNGSCHR_inferChrSize, NOSITE){
       StandardChroms ourChroms;
       ourChroms.emptyChr.inferChrSize();
       EXPECT_EQ(0,ourChroms.emptyChr.getChromSize() );
 }
 TEST(uBasicNGSCHR_inferChrSize, MANYSITE){
       StandardChroms ourChroms;
       ourChroms.manyChr.inferChrSize();
       EXPECT_EQ(300,ourChroms.manyChr.getChromSize() );
 }
/**<  Count Unique*/
 TEST(uBasicNGSCHR_countUnique, ONESITE){
       StandardChroms ourChroms;
       EXPECT_EQ(1,ourChroms.oneChr.countUnique() );
 }
TEST(uBasicNGSCHR_countUnique, NOSITE){
         StandardChroms ourChroms;
         EXPECT_EQ(0 , ourChroms.emptyChr.countUnique());
 }
 TEST(uBasicNGSCHR_countUnique, MANYSITENOUNIQUE){
       StandardChroms ourChroms;
       EXPECT_EQ( 3, ourChroms.manyChr.countUnique());
 }
TEST(uBasicNGSCHR_countUnique, MANYSITEWITHUNIQUE){
       uBasicNGSChrom newManyChr("chr1");
       newManyChr.addData(uBasicNGS("chr1",100,200));
       newManyChr.addData(uBasicNGS("chr1",100,200));
       newManyChr.addData(uBasicNGS("chr1",150,200));
       newManyChr.addData(uBasicNGS("chr1",100,2500));
       newManyChr.addData(uBasicNGS("chr1",230,300));
       newManyChr.addData(uBasicNGS("chr1",120,250));
       EXPECT_EQ(5 , newManyChr.countUnique());
 }

/**< Testing DivideItemsIntoNBin */
// TODO: strict
// TODO: throw test
TEST_F(ChromDivide, DIVIDEINTONBINFAILCHROM){
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(7));
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(32));
}
TEST_F(ChromDivide, DIVIDEINTONBINCHROMIGNORE){
    //std::cout << uChromTestOverlap.count();
    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::IGNORE);
    EXPECT_EQ(uChromTestOverlap.count(),12);
}
TEST_F(ChromDivide, DIVIDEINTOBINEXTEND){

    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::EXTEND);
    EXPECT_EQ(uChromTestOverlap.count(),12);
}
TEST_F(ChromDivide, DIVIDEINTOBINADD){

    uChromTestOverlap.divideItemsIntoNBins(4, SplitType::ADD);
    EXPECT_EQ(uChromTestOverlap.count(),13);
}
TEST_F(ChromDivide, DIVIDECHROMINTONBINOFSIZEFAIL){
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoBinofSize(30));
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoBinofSize(300));
}

// TODO: strict
// TODO: throw test
TEST_F(ChromDivide, DIVIDECHROMINTONBINOFSIZEIGNORE){

    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::IGNORE);
    EXPECT_EQ(uChromTestOverlap.count(),5);
}
TEST_F(ChromDivide, DIVIDECHROMINTOBINOFSIZEEXTEND){
    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::EXTEND);
    EXPECT_EQ(uChromTestOverlap.count(),5);
}
TEST_F(ChromDivide, DIVIDECHROMINTOBINOFSIZEADD){
    uChromTestOverlap.divideItemsIntoBinofSize(90, SplitType::ADD);

    EXPECT_EQ(uChromTestOverlap.count(),8);
}

/**<  Random generation test*/
TEST(uBasicNGSCHR_generateRandomSite, NOCHRSIZE){
       StandardChroms ourChroms;
       std::random_device rd;
       std::mt19937 gen(rd());
       EXPECT_THROW(ourChroms.manyChr.generateRandomSite(100,gen,0,"test"),param_throw );
 }
TEST(uBasicNGSCHR_generateRandomSite, LARGETTHENCHRSIZE){
       StandardChroms ourChroms;
       ourChroms.manyChr.setChromSize(500);
       std::random_device rd;
       std::mt19937 gen(rd());
       EXPECT_THROW(ourChroms.manyChr.generateRandomSite(1000,gen,0,"test"),param_throw );
 }

TEST(uBasicNGSCHR_generateRandomSite, MAKEITEM){
       StandardChroms ourChroms;
       ourChroms.manyChr.setChromSize(500);
       std::random_device rd;
       std::mt19937 gen(rd());
       uBasicNGS aItem= ourChroms.manyChr.generateRandomSite(100,gen,0);
       EXPECT_EQ(100,aItem.getLength());
       uBasicNGS aItem2= ourChroms.manyChr.generateRandomSite(150,gen,0);
       EXPECT_EQ(150,aItem2.getLength());
 }
TEST(uBasicNGSCHR_generateRandomSite, MAKEMANYDIFFERENTITEMS){
       StandardChroms ourChroms;
       ourChroms.manyChr.setChromSize(500);

       std::mt19937 gen(432);
       uBasicNGS aItem= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uBasicNGS aItem2= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uBasicNGS aItem3= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       uBasicNGS aItem4= ourChroms.manyChr.generateRandomSite(100,gen,0,"test");
       EXPECT_EQ(aItem2.getLength(),aItem.getLength());
       EXPECT_EQ(aItem3.getLength(),aItem2.getLength());
       EXPECT_EQ(aItem4.getLength(),aItem3.getLength());

       EXPECT_NE(aItem.getStart(),aItem2.getStart());
       EXPECT_NE(aItem2.getStart(),aItem3.getStart());
       EXPECT_NE(aItem3.getStart(),aItem4.getStart());
 }
TEST(uBasicNGSCHR_generateRandomSite, MAKEMANYEXCLUSIONLIST){
      StandardChroms ourChroms;
      uBasicNGSChrom exclusionChrom("chr1");
      exclusionChrom.addData(uBasicNGS("chr1",100,200));
      exclusionChrom.addData(uBasicNGS("chr1",150,500));
      exclusionChrom.addData(uBasicNGS("chr1",800,1000));
      exclusionChrom.setChromSize(1500);
      ourChroms.emptyChr.setChromSize(1500);
      std::mt19937 gen(432);
      for(size_t i =0; i<=2000; i++)
        ourChroms.emptyChr.addData(ourChroms.emptyChr.generateRandomSite(100,gen,exclusionChrom,4,"test"));


      auto chromOverlapping= ourChroms.emptyChr.getOverlapping(exclusionChrom);
      EXPECT_EQ(chromOverlapping.count(),0);
 }

TEST(uBasicNGSCHR_addNRandomSite, MAKECORRECTNUMBER){
       StandardChroms ourChroms;
       ourChroms.emptyChr.setChromSize(1000);
       std::random_device rd;
       std::mt19937 gen(rd());
       ourChroms.emptyChr.addNRandomSite(50,10,gen,0);
       EXPECT_EQ(10,ourChroms.emptyChr.count());
 }
 TEST(uBasicNGSCHR_addNRandomSite, MAKECORRECTNUMBERMANYDIFFERENT){
       StandardChroms ourChroms;
       ourChroms.emptyChr.setChromSize(1000);
       std::random_device rd;
       std::mt19937 gen(4523);
       ourChroms.emptyChr.addNRandomSite(50,10,gen,0);
       EXPECT_EQ(10,ourChroms.emptyChr.count());
 }
 TEST(uBasicNGSCHR_addNRandomSite, MAKEMANYEXCLUSIONLIST){
   StandardChroms ourChroms;
      uBasicNGSChrom exclusionChrom("chr1");
      exclusionChrom.addData(uBasicNGS("chr1",100,200));
      exclusionChrom.addData(uBasicNGS("chr1",150,500));
      exclusionChrom.addData(uBasicNGS("chr1",800,1000));
      exclusionChrom.setChromSize(1500);
      ourChroms.emptyChr.setChromSize(1500);
      std::mt19937 gen(432);
      ourChroms.emptyChr.addNRandomSite(100,2000,gen,exclusionChrom,4);
      EXPECT_EQ(ourChroms.emptyChr.count(),2000);
      auto chromOverlapping= ourChroms.emptyChr.getOverlapping(exclusionChrom);
      EXPECT_EQ(chromOverlapping.count(),0);
 }

/**< Get overlapping */
 TEST(uBasicNGSCHR_getOverlapping, NORMAL){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",50,100));
       uBasicNGSChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(1,retChr.count());
 }
 TEST(uBasicNGSCHR_getOverlapping, FIRSTEMPTY){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",50,100));
       uBasicNGSChrom retChr=ourChroms.emptyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uBasicNGSCHR_getOverlapping, SECONDEMPTY){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       uBasicNGSChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uBasicNGSCHR_getOverlapping, SOMEOVERLAP){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",0,220));
       uBasicNGSChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(2,retChr.count());
 }
  TEST(uBasicNGSCHR_getOverlapping, ALLOVERLAP){
       StandardChroms ourChroms;
       uBasicNGSChrom twoItems("chr1");
       twoItems.addData(uBasicNGS("chr1",0,220));
       twoItems.addData(uBasicNGS("chr1",230,420));
       uBasicNGSChrom retChr=ourChroms.manyChr.getOverlapping(twoItems);
       EXPECT_EQ(3,retChr.count());
 }
 TEST(uBasicNGSCHR_getOverlapping, NOOVERLAP){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",500,1000));
       uBasicNGSChrom retChr=ourChroms.emptyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
 TEST(uBasicNGSCHR_getOverlapping, DIFFERENTCHR){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr2");
       oneItem.addData(uBasicNGS("chr2",0,1000));
       uBasicNGSChrom retChr=ourChroms.manyChr.getOverlapping(oneItem);
       EXPECT_EQ(0,retChr.count());
 }
/**< Get overlapping Count */

 TEST(uBasicNGSCHR_getOverlappingCount, NORMAL){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",50,100));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),1);
 }
 TEST(uBasicNGSCHR_getOverlappingCount, FIRSTEMPTY){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",50,100));
       EXPECT_EQ(ourChroms.emptyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uBasicNGSCHR_getOverlappingCount, SECONDEMPTY){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uBasicNGSCHR_getOverlappingCount, SOMEOVERLAP){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",0,220));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),2);
 }
  TEST(uBasicNGSCHR_getOverlappingCount, ALLOVERLAP){
       StandardChroms ourChroms;
       uBasicNGSChrom twoItems("chr1");
       twoItems.addData(uBasicNGS("chr1",0,220));
       twoItems.addData(uBasicNGS("chr1",230,420));
       EXPECT_EQ(3,ourChroms.manyChr.getOverlappingCount(twoItems));


 }
 TEST(uBasicNGSCHR_getOverlappingCount, NOOVERLAP){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",500,1000));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }
 TEST(uBasicNGSCHR_getOverlappingCount, DIFFERENTCHR){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr2");
       oneItem.addData(uBasicNGS("chr2",0,1000));
       EXPECT_EQ(ourChroms.manyChr.getOverlappingCount(oneItem),0);
 }

/**< Get not Overlapping */
 TEST(uBasicNGSCHR_getNotOverlapping, NORMAL){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",50,100));
       uBasicNGSChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),2);
 }
 TEST(uBasicNGSCHR_getNotOverlapping, FIRSTEMPTY){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",50,100));
       uBasicNGSChrom retChr=ourChroms.emptyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),0);
 }
  TEST(uBasicNGSCHR_getNotOverlapping, SECONDEMPTY){
       StandardChroms ourChroms;
       uBasicNGSChrom noItems("chr1");
       uBasicNGSChrom retChr=ourChroms.manyChr.getNotOverlapping(noItems);
       EXPECT_EQ(retChr.count(),3);
 }
  TEST(uBasicNGSCHR_getNotOverlapping, SOMEOVERLAP){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",0,220));
       uBasicNGSChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),1);
 }
  TEST(uBasicNGSCHR_getNotOverlapping, NOOVERLAP){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr1");
       oneItem.addData(uBasicNGS("chr1",500,1000));
       uBasicNGSChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),3);
 }
 TEST(uBasicNGSCHR_getNotOverlapping, DIFFERENTCHR){
       StandardChroms ourChroms;
       uBasicNGSChrom oneItem("chr2");
       oneItem.addData(uBasicNGS("chr2",0,1000));
       uBasicNGSChrom retChr=ourChroms.manyChr.getNotOverlapping(oneItem);
       EXPECT_EQ(retChr.count(),3);
 }

/**<  getDistinct && getDistinctCount &&removeDistinct*/
TEST(uBasicNGSCHR_getDistinct, EMPTY)
{
     uBasicNGSChrom emptyChr("chr2");
     uBasicNGSChrom testChrom =emptyChr.getDistinct(500, 2000);
     EXPECT_EQ(0,testChrom.count());
     EXPECT_EQ(0,emptyChr.getDistinctCount(500, 2000));
}
TEST(uBasicNGSCHR_getDistinct, MULTIPLE)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     uBasicNGSChrom testChrom =ourChroms.manyChr.getDistinct(500, 2000);
     EXPECT_EQ(3,testChrom.count());
     EXPECT_EQ(3,ourChroms.manyChr.getDistinctCount(500, 2000));
}
TEST(uBasicNGSCHR_getDistinct, GETNONE)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     uBasicNGSChrom testChrom =ourChroms.manyChr.getDistinct(0, 2000);
     EXPECT_EQ(0,testChrom.count());
     EXPECT_EQ(0,ourChroms.manyChr.getDistinctCount(0, 2000));
}
TEST(uBasicNGSCHR_getDistinct, CUSTOMSORT)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     uBasicNGSChrom testChrom =ourChroms.manyChr.getDistinct(80, 130);
     EXPECT_EQ(2,testChrom.count());
     EXPECT_EQ(2,ourChroms.manyChr.getDistinctCount(80, 130));
}
TEST(uBasicNGSCHR_getDistinct, FAILNOSORT)
{
     StandardChroms ourChroms;
     EXPECT_THROW(ourChroms.manyChr.getDistinct(80, 130),unsorted_throw );
     EXPECT_THROW(ourChroms.manyChr.getDistinctCount(80, 130),unsorted_throw );
}
TEST(uBasicNGSCHR_getDistinct, MORECUSTOMSORT)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);

     uBasicNGSChrom testChrom =ourChroms.manyChr.getDistinct(80, 130);
     EXPECT_EQ(2,testChrom.count());
     EXPECT_EQ(2,ourChroms.manyChr.getDistinctCount(80, 130));
}


/**< RemoveDistinct */
TEST(uBasicNGSCHR_removeDistinct, EMPTY)
{
     uBasicNGSChrom emptyChr("chr2");
     uBasicNGSChrom testChrom =emptyChr.removeDistinct(500, 2000);
     EXPECT_EQ(0,testChrom.count());
     EXPECT_EQ(0,emptyChr.count());
}
TEST(uBasicNGSCHR_removeDistinct, MULTIPLE)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     uBasicNGSChrom testChrom =ourChroms.manyChr.removeDistinct(500, 2000);
     EXPECT_EQ(3,testChrom.count());
     EXPECT_EQ(0,ourChroms.manyChr.count());
}
TEST(uBasicNGSCHR_removeDistinct, GETNONE)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     uBasicNGSChrom testChrom =ourChroms.manyChr.removeDistinct(0, 2000);
     EXPECT_EQ(0,testChrom.count());
     EXPECT_EQ(3,ourChroms.manyChr.count());
}
TEST(uBasicNGSCHR_removeDistinct, CUSTOMSORT)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     uBasicNGSChrom testChrom =ourChroms.manyChr.removeDistinct(80, 130);
     EXPECT_EQ(2,testChrom.count());
     EXPECT_EQ(1,ourChroms.manyChr.count());
}
TEST(uBasicNGSCHR_removeDistinct, FAILNOSORT)
{
     StandardChroms ourChroms;
     EXPECT_THROW(ourChroms.manyChr.removeDistinct(80, 130),unsorted_throw );
}

/**< Test getSubset, RemoveSubset */
TEST(uBasicNGSCHR_getSubset, EMPTY)
{
     StandardChroms ourChroms;
     uBasicNGSChrom testChrom =ourChroms.emptyChr.getSubset(0, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uBasicNGSCHR_getSubset, MANY)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     uBasicNGSChrom testChrom =ourChroms.manyChr.getSubset(0, 2000);
     EXPECT_EQ(3,testChrom.count());
}
TEST(uBasicNGSCHR_getSubset, UNSORTED)
{
     StandardChroms ourChroms;
     EXPECT_THROW(ourChroms.manyChr.getSubset(80, 130),unsorted_throw );
}
TEST(uBasicNGSCHR_getSubset, CUSTOMEMPTY)
{
     StandardChroms ourChroms;
     ourChroms.emptyChr.sortSites(ourChroms.emptyChr.compareLength,&uBasicNGS::getLength);
     uBasicNGSChrom testChrom =ourChroms.emptyChr.getSubset(80, 130);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uBasicNGSCHR_getSubset, CUSTOMESOME)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     uBasicNGSChrom testChrom =ourChroms.manyChr.getSubset(80, 130);
     EXPECT_EQ(1,testChrom.count());
}

TEST(uBasicNGSCHR_removeSubset, CUSTOMBORDER)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     uBasicNGSChrom testChrom =ourChroms.manyChr.removeSubset(80, 130);

     EXPECT_EQ(2,ourChroms.manyChr.count());
     EXPECT_EQ(1,testChrom.count());
}

TEST(uBasicNGSCHR_removeSubset, CUSTOMALL)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     uBasicNGSChrom testChrom =ourChroms.manyChr.removeSubset(0, 200);
     EXPECT_EQ(3,testChrom.count());
}
TEST(uBasicNGSCHR_removeSubset, EMPTY)
{
     StandardChroms ourChroms;
     uBasicNGSChrom testChrom =ourChroms.emptyChr.removeSubset(0, 2000);
     EXPECT_EQ(0,testChrom.count());
}
TEST(uBasicNGSCHR_removeSubset, ALLREMOVED)
{

     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     uBasicNGSChrom testChrom =ourChroms.manyChr.removeSubset(0, 2000);
     EXPECT_EQ(0,ourChroms.manyChr.count());
     EXPECT_EQ(3,testChrom.count());
}
TEST(uBasicNGSCHR_removeSubset, SOMEREMOVED)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     uBasicNGSChrom testChrom =ourChroms.manyChr.removeSubset(0, 150);
     EXPECT_EQ(1,ourChroms.manyChr.count());
     EXPECT_EQ(2,testChrom.count());

}

/**<  getSubsetCount*/

TEST(uBasicNGSCHR_getSubsetCount, NORMAL)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     EXPECT_EQ(3,ourChroms.manyChr.getSubsetCount(0, 2000));
}

TEST(uBasicNGSCHR_getSubsetCount, MIDDLE)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites();
     EXPECT_EQ(2,ourChroms.manyChr.getSubsetCount(125, 300));
}

TEST(uBasicNGSCHR_getSubsetCount, EMPTY)
{
     StandardChroms ourChroms;
     EXPECT_EQ(0,ourChroms.emptyChr.getSubsetCount(0, 2000));
}

TEST(uBasicNGSCHR_getSubsetCount, MANYCUSTOMSORT)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     EXPECT_EQ(2,ourChroms.manyChr.getSubsetCount(40, 120));
}
/**< AddData */
// TODO: With token
TEST(uBasicNGSCHR_addData, VALID)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     ourChroms.manyChr.addData(uBasicNGS("chr1", 2032,3245));
     EXPECT_FALSE(ourChroms.manyChr.getSortedStatus());
     EXPECT_EQ(4,ourChroms.manyChr.count());

}
TEST(uBasicNGSCHR_addData, INVALIDCHR)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     EXPECT_THROW(ourChroms.manyChr.addData(uBasicNGS("chr2", 2032,3245)),ugene_exception_base);
}


/**<  GetStartFunct, getEndFunct */
TEST(uBasicNGSCHR_getStartFunct, HASBEENSET)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
     std::function<float(const uBasicNGS*)>  my_funct= &uBasicNGS::getLength;
     EXPECT_NE(nullptr,ourChroms.manyChr.getStartFunct());
}
TEST(uBasicNGSCHR_getStartFunct, NOTBEENSET)
{
     StandardChroms ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getStartFunct());
}

TEST(uBasicNGSCHR_getEndFunct, HASBEENSET)
{
     StandardChroms ourChroms;
     ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength,&uBasicNGS::getLength);
     std::function<float(const uBasicNGS*)>  my_funct= &uBasicNGS::getLength;
     EXPECT_NE(nullptr,ourChroms.manyChr.getEndFunct());
}

TEST(uBasicNGSCHR_getEndFunct, NOTBEENSET)
{
     StandardChroms ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getEndFunct());
}

TEST(uBasicNGSCHR_getCompFunct, BEENSET)
{
      StandardChroms ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength,&uBasicNGS::getLength);
      EXPECT_NE(nullptr,ourChroms.manyChr.getCompFunct());
}
TEST(uBasicNGSCHR_getCompFunct, NOTBEENSET)
{
   StandardChroms ourChroms;
     EXPECT_NE(nullptr,ourChroms.manyChr.getCompFunct());
}

/**<  setChromSize*/

TEST(uBasicNGSCHR_setChromSize, VALID)
{
     StandardChroms ourChroms;
     EXPECT_NO_THROW(ourChroms.emptyChr.setChromSize(20000));
     EXPECT_EQ(20000,ourChroms.emptyChr.getChromSize());
}
TEST(uBasicNGSCHR_setChromSize, INVALID_UNDER0)
{
     StandardChroms ourChroms;
     EXPECT_THROW(ourChroms.emptyChr.setChromSize(-200),param_throw);
}

/**< getSite */
TEST(uBasicNGSCHR_getSite, VALIDREQUEST)
{
     StandardChroms ourChroms;
     uBasicNGS oneItem;
     EXPECT_NO_THROW(oneItem=ourChroms.manyChr.getSite(2));
     EXPECT_TRUE(oneItem.isEqual(uBasicNGS("chr1",120,250)));

}
TEST(uBasicNGSCHR_getSite, INVALID)
{
     StandardChroms ourChroms;
     uBasicNGS oneItem;
     EXPECT_THROW(oneItem=ourChroms.manyChr.getSite(20),param_throw);
}

/**< SortSites */
TEST(uBasicNGSCHR_sortSites, DEFAULT)
{
      StandardChroms ourChroms;
      ourChroms.manyChr.sortSites();
      auto itr= ourChroms.manyChr.begin();
      //Check manually the order
      EXPECT_TRUE((itr)->isEqual(uBasicNGS("chr1",100,200)));
      EXPECT_TRUE((itr+1)->isEqual(uBasicNGS("chr1",120,250)));
      EXPECT_TRUE((itr+2)->isEqual(uBasicNGS("chr1",230,300)));
}

TEST(uBasicNGSCHR_sortSites, CUSTOM)
{
      StandardChroms ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength);
      //Check manually the order
      auto itr= ourChroms.manyChr.begin();

      EXPECT_TRUE((itr)->isEqual(uBasicNGS("chr1",230,300)));
      EXPECT_TRUE((itr+1)->isEqual(uBasicNGS("chr1",100,200)));
      EXPECT_TRUE((itr+2)->isEqual(uBasicNGS("chr1",120,250)));
}

// TODO: Test throw when not sorted
// TODO: Test after last element or before first element
TEST(uBasicNGSCHR_findNext, STANDARD)
{
      StandardChroms ourChroms;
      ourChroms.manyChr.sortSites();
      auto first=ourChroms.manyChr.findNextSite(195);
      EXPECT_TRUE(first->isEqual(uBasicNGS("chr1",230,300)));
}
TEST(uBasicNGSCHR_findNext, CUSTOM)
{
      StandardChroms ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength,&uBasicNGS::getLength);
      auto first=ourChroms.manyChr.findNextSite(125);
      EXPECT_TRUE(first->isEqual(uBasicNGS("chr1",120,250)));
      first=ourChroms.manyChr.findNextSite(0);
      EXPECT_TRUE(first->isEqual(uBasicNGS("chr1",230,300)));
}

TEST(uBasicNGSCHR_findNext, EMPTY)
{
      StandardChroms ourChroms;
      ourChroms.emptyChr.sortSites();
      auto first=ourChroms.emptyChr.findNextSite(120);
      EXPECT_EQ(first,ourChroms.emptyChr.end());
}

TEST(uBasicNGSCHR_findPrec, STANDARD)
{
      StandardChroms ourChroms;
      ourChroms.manyChr.sortSites();
      auto first=ourChroms.manyChr.findPrecedingSite(195);
      EXPECT_TRUE(first->isEqual(uBasicNGS("chr1",120,250)));
}

TEST(uBasicNGSCHR_findPrec, CUSTOM)
{
      StandardChroms ourChroms;
      ourChroms.manyChr.sortSites(ourChroms.manyChr.compareLength,&uBasicNGS::getLength,&uBasicNGS::getLength);
      auto first=ourChroms.manyChr.findPrecedingSite(125);
      EXPECT_TRUE(first->isEqual(uBasicNGS("chr1",100,200)));
      first=ourChroms.manyChr.findPrecedingSite(400);
      EXPECT_TRUE(first->isEqual(uBasicNGS("chr1",120,250)));
}


TEST(uBasicNGSCHR_findPrec, EMPTY)
{
      StandardChroms ourChroms;
      ourChroms.emptyChr.sortSites();
      auto first=ourChroms.emptyChr.findPrecedingSite(120);
      EXPECT_EQ(first,ourChroms.emptyChr.end());
}
/**< isSorted() */
// TODO: Check if function that possibly modify sorted status does so correctly
TEST(uBasicNGSCHR_isSorted, NOT)
{
     StandardChroms ourChroms;
     EXPECT_FALSE(ourChroms.manyChr.isSorted());
}
TEST(uBasicNGSCHR_isSorted, EMPTY)
{
       StandardChroms ourChroms;
       EXPECT_TRUE(ourChroms.emptyChr.isSorted());
}
TEST(uBasicNGSCHR_isSorted, NORMAL)
{
       StandardChroms ourChroms;
       ourChroms.manyChr.sortSites();
       EXPECT_TRUE(ourChroms.emptyChr.isSorted());
       EXPECT_FALSE(ourChroms.manyChr.isSorted(uBasicNGSChrom::compareLength));
}
TEST(uBasicNGSCHR_isSorted, CUSTOM)
{
       StandardChroms ourChroms;
       ourChroms.manyChr.sortSites(uBasicNGSChrom::compareLength);
       EXPECT_TRUE(ourChroms.manyChr.isSorted());
       EXPECT_FALSE(ourChroms.manyChr.isSorted(uBasicNGSChrom::compareStart));
}

TEST(uBasicNGSCHR_RemoveSite, SINGLEITERATOR)
{
       StandardChroms ourChroms;
       ourChroms.manyChr.removeSite(ourChroms.manyChr.begin());
       EXPECT_EQ(2,ourChroms.manyChr.count());
       EXPECT_TRUE(ourChroms.manyChr.begin()->isEqual(uBasicNGS("chr1",230,300)));
        StandardChroms newChroms;
        auto itr =newChroms.manyChr.begin();
        itr++;
        newChroms.manyChr.removeSite(itr);
        EXPECT_EQ(2,newChroms.manyChr.count());
        EXPECT_TRUE(newChroms.manyChr.begin()->isEqual(uBasicNGS("chr1",100,200)));
}
TEST(uBasicNGSCHR_RemoveSite, ITR_RANGE)
{
       StandardChroms ourChroms;
       ourChroms.manyChr.removeSite(ourChroms.manyChr.begin(), ourChroms.manyChr.end());
       EXPECT_EQ(0,ourChroms.manyChr.count());
        StandardChroms newChroms;
        auto itr =newChroms.manyChr.begin();
        auto itr2 =newChroms.manyChr.begin();
        itr2+=2;
        newChroms.manyChr.removeSite(itr,itr2);
        EXPECT_EQ(1,newChroms.manyChr.count());
        EXPECT_TRUE(newChroms.manyChr.begin()->isEqual(uBasicNGS("chr1",120,250)));
}
/**< writeWithWriter */

TEST(uBasicNGSCHR_writeWithWriter, NORMALBED)
{
       StandardChroms ourChroms;
       uWriter bedWriter(&cout,"BED4");
       EXPECT_NO_THROW(ourChroms.manyChr.writeWithWriter(bedWriter));
}
