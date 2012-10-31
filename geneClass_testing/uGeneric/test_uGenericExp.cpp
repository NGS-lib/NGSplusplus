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

class GeneExp : public uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>, uGenericNGS> {

public:
  GeneExp():uGenericNGSExperiment(){};
  //ourDerivedClass(std::string chrom):uGenericNGSExperiment(chrom){};
  //ourDerivedClass(std::string chrom, long int size):uGenericNGSExperiment(chrom,size){};
};


class TestExp : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {
    uExpTest.addSite(uGenericNGS("chr1", 300, 500));
    uExpTest.addSite(uGenericNGS("chr1", 100, 199));
    uExpTest.addSite(uGenericNGS("chr1", 100, 299));
  }
GeneExp uExpTest;
};

/**< Start chrom test */
TEST_F(TestExp, SORTEXP){
    /**< Unsorted, fail */
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites();
    EXPECT_TRUE(uExpTest.isSorted());
    uExpTest.addSite(uGenericNGS("chr1", 2, 50));
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites();
    EXPECT_TRUE(uExpTest.isSorted());
}

/**< Start chrom test */
TEST_F(TestExp, SORTCUSTOM){
    /**< Unsorted, fail */
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites();
    EXPECT_TRUE(uExpTest.isSorted());
    uExpTest.addSite(uGenericNGS("chr1", 2, 50));
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites();
    EXPECT_TRUE(uExpTest.isSorted());
}
