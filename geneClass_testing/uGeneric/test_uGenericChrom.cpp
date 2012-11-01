
#include "uFormats.h"
#include "uTags.h"
#include "uRegion.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include "utility/utility.h"
#include <time.h>
#include "gtest.h"
using namespace NGS;

class ourDerivedClass : public uGenericNGSChrom<uGenericNGS> {

public:
  ourDerivedClass():uGenericNGSChrom(){};
  ourDerivedClass(std::string chrom):uGenericNGSChrom(chrom){};
  ourDerivedClass(std::string chrom, long int size):uGenericNGSChrom(chrom,size){};
};

class ourDerivedExperiment : public uGenericNGSExperiment<ourDerivedClass, uGenericNGS> {

};

// To use a test fixture, derive a class from testing::Test.
class ChromDivide : public testing::Test {

 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {
    uChromTestOverlap.addSite(uGenericNGS("chr1", 300, 500));
    uChromTestOverlap.addSite(uGenericNGS("chr1", 100, 199));
    uChromTestOverlap.addSite(uGenericNGS("chr1", 100, 299));

  }
ourDerivedClass uChromTestOverlap;
};

/**< Start chrom test */
TEST_F(ChromDivide, SORT){

    /**< Unsorted, fail */
    EXPECT_FALSE(uChromTestOverlap.isSorted());
    uChromTestOverlap.sortSites();
    EXPECT_TRUE(uChromTestOverlap.isSorted());
    uChromTestOverlap.addSite(uGenericNGS("chr1", 200, 800));
    uChromTestOverlap.addSite(uGenericNGS("chr1", 100, 800));
    uChromTestOverlap.addSite(uGenericNGS("chr1", 300, 800));
    EXPECT_FALSE(uChromTestOverlap.isSorted());
}


TEST_F(ChromDivide, FIND_TEST){

    /**< Unsorted, fail */

    uChromTestOverlap.addSite(uGenericNGS("chr1", 200, 800));
    uChromTestOverlap.addSite(uGenericNGS("chr1", 250, 800));
     uChromTestOverlap.sortSites();
    auto first=uChromTestOverlap.findPrecedingSite(195);
    EXPECT_EQ(100,first->getStart());

    auto second=uChromTestOverlap.findNextSite(235);
    EXPECT_EQ(250,second->getStart());

    EXPECT_EQ(uChromTestOverlap.end(), uChromTestOverlap.findNextSite(1500));
   //=uChromTestOverlap.findPrecedingSite(195);

}

/**< Testing DivideItemsIntoNBin */
TEST_F(ChromDivide, DIVIDEINTONBINFAILCHROM){
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(7));
    EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(32));
}

TEST_F(ChromDivide, DIVIDEINTONBINCHROMIGNORE){

  //  cout << uChromTestOverlap.count();
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

  //  cout << uChromTestOverlap.count();
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
