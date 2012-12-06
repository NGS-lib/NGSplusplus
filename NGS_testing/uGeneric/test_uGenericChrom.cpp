
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
#include "gtest.h"
using namespace NGS;

//
//class uBasicNGS:public uGenericNGS<uBasicNGS>{
//
//    public:
//    uBasicNGS():uGenericNGS(){};
//    uBasicNGS(std::string chr, int start,int end):uGenericNGS(chr,start,end){};
//    uBasicNGS(uGenericNGS otherItem):uGenericNGS(otherItem)
//    { };
//    uBasicNGS& operator=(const uBasicNGS& copFrom)=default;
//    uBasicNGS(const uBasicNGS&)=default;
//};
//
//class ourDerivedClass : public uGenericNGSChrom<ourDerivedClass,uBasicNGS> {
//
//public:
//  ourDerivedClass():uGenericNGSChrom(){};
//  ourDerivedClass(std::string chrom):uGenericNGSChrom(chrom){};
//  ourDerivedClass(std::string chrom, long int size):uGenericNGSChrom(chrom,size){};
//
//  ourDerivedClass& operator=(const ourDerivedClass& copFrom)=default;
//  ourDerivedClass(const ourDerivedClass&)=default;
//
//};
//
//class ourDerivedExperiment : public uGenericNGSExperiment<ourDerivedExperiment,ourDerivedClass, uBasicNGS> {
//
//};
// To use a test fixture, derive a class from testing::Test.
class ChromDivide : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {

    uChromTestOverlap.setChr("chr1");
    uChromTestOverlap.addSite(uBasicNGS ("chr1", 300, 500));
    uChromTestOverlap.addSite(uBasicNGS("chr1", 100, 199));
    uChromTestOverlap.addSite(uBasicNGS("chr1", 100, 299));
  }
uBasicNGSChrom uChromTestOverlap;
};

class ChromSubset : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {
    uChromSubsetTest.setChr("chr1");
    uChromSubsetTest.addSite(uBasicNGS("chr1", 300, 500));
    uChromSubsetTest.addSite(uBasicNGS("chr1", 100, 199));
    uChromSubsetTest.addSite(uBasicNGS("chr1", 100, 299));
    uChromSubsetTest.sortSites();
  }
uBasicNGSChrom uChromSubsetTest;
};


/**< Start chrom test */
TEST_F(ChromDivide, SORT){

    /**< Unsorted, fail */
    EXPECT_FALSE(uChromTestOverlap.isSorted());
    uChromTestOverlap.sortSites();
    EXPECT_TRUE(uChromTestOverlap.isSorted());
    uChromTestOverlap.addSite(uBasicNGS("chr1", 200, 800));
    uChromTestOverlap.addSite(uBasicNGS("chr1", 100, 800));
    uChromTestOverlap.addSite(uBasicNGS("chr1", 300, 800));
    EXPECT_FALSE(uChromTestOverlap.isSorted());
}


TEST_F(ChromDivide, FIND_TEST){
    uChromTestOverlap.addSite(uBasicNGS("chr1", 200, 800));
    uChromTestOverlap.addSite(uBasicNGS("chr1", 250, 800));
    uChromTestOverlap.sortSites();
    auto first=uChromTestOverlap.findPrecedingSite(195);
    EXPECT_EQ(100,first->getStart());

    auto second=uChromTestOverlap.findNextSite(235);
    EXPECT_EQ(250,second->getStart());

    EXPECT_EQ(uChromTestOverlap.end(), uChromTestOverlap.findNextSite(1500));

}

/**< Testing DivideItemsIntoNBin */
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
/**< Testing DivideItemsIntoBinofSize */
TEST_F(ChromDivide, DIVIDECHROMINTONBINOFSIZEFAIL){
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoBinofSize(30));

    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoBinofSize(300));
}

TEST_F(ChromDivide, DIVIDECHROMINTONBINOFSIZEIGNORE){

//    cout << uChromTestOverlap.count();
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


TEST_F(ChromSubset, distinctChrTestBorder)
{
     auto testChrom=uChromSubsetTest.getDistinct(5, 200);
     EXPECT_EQ(1,testChrom.count());
}

TEST_F(ChromSubset, subsetCountChr)
{

     EXPECT_EQ(3,uChromSubsetTest.getSubsetCount(1, 2000));
}

TEST(ChromSubsets, distinctChrTestEmpty)
{
    //uTagsChrom emptyChr("chr2");
   // auto testChrom =emptyChr.getDistinct(500, 2000);
     uBasicNGSChrom emptyChr("chr2");
     auto testChrom =emptyChr.getDistinct(500, 2000);
     EXPECT_EQ(0,testChrom.count());
}




