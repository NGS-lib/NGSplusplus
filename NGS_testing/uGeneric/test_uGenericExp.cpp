#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include "NGS++.h"
#include <time.h>
#include <functional>
#include "gtest.h"
using namespace NGS;
using namespace std;
class GeneExp : public uGenericNGSExperiment<GeneExp,uBasicNGSChrom, uBasicNGS> {


public:
  GeneExp():uGenericNGSExperiment(){};
  //ourDerivedClass(std::string chrom):uGenericNGSExperiment(chrom){};
  //ourDerivedClass(std::string chrom, long int size):uGenericNGSExperiment(chrom,size){};
};

class TestExp : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {
    uExpTest.addSite(uGenericNGS("chr1", 300, 500, 0.1f));
    uExpTest.addSite(uGenericNGS("chr1", 100, 199,0.5f));
    uExpTest.addSite(uGenericNGS("chr1", 100, 299,0.3f));
  }
GeneExp uExpTest;
};

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


TEST_F(TestExp, SORTCUSTOM){

    auto comp=[](const uGenericNGS & item1, const uGenericNGS & item2){
             return item1.getScore() < item2.getScore();
    };
    /**< Unsorted, fail */
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites(comp);
    EXPECT_TRUE(uExpTest.isSorted());
    uExpTest.addSite(uGenericNGS("chr1", 2, 50, 0.05f));
    uExpTest.addSite(uGenericNGS("chr1", 2, 50, 0.7f));
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites();
    EXPECT_TRUE(uExpTest.isSorted());
}


TEST_F(TestExp, GETSUBSET){

    uExpTest.addSite(uGenericNGS("chr1", 2, 50, 0.05f));
    uExpTest.addSite(uGenericNGS("chr1", 5, 25, 0.7f));
    uExpTest.sortSites();
    auto resultChrom=uExpTest.getSubset("chr1",99,200);
    EXPECT_EQ(resultChrom.count(), 2);

}


TEST_F(TestExp, GETSUBSETLARGE){

    uExpTest.addSite(uGenericNGS("chr1", 180000050, 180500000 ));
    uExpTest.addSite(uGenericNGS("chr1", 180050000, 180500000));
    uExpTest.addSite(uGenericNGS("chr1", 180550000, 180900000));
    uExpTest.addSite(uGenericNGS("chr1", 5, 25, 0.7f));
    uExpTest.addSite(uGenericNGS("chr1", 5, 25, 0.7f));

    uExpTest.sortSites();
    auto resultChrom=uExpTest.getSubset("chr1",180000000,181000000);
    EXPECT_EQ(resultChrom.count(), 3);

}


TEST_F(TestExp, GETSUBSETCUSTOM){

  auto comp=[](const uGenericNGS & item1, const uGenericNGS & item2){
             return item1.getScore() < item2.getScore();
    };

    uExpTest.addSite(uGenericNGS("chr1", 2, 50, 0.05f));
    uExpTest.addSite(uGenericNGS("chr1", 5, 25, 0.7f));
    auto getScoreFunct=std::mem_fn((float(uGenericNGS::*)()const)&uGenericNGS::getScore);

    uExpTest.sortSites(comp,getScoreFunct,getScoreFunct);
    auto resultChrom=uExpTest.getSubset("chr1",0.4f,0.9f);
    EXPECT_EQ(resultChrom.count(), 3);
}


TEST(TestExpNOCASE, LOADBED){

   // uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>, uGenericNGS>  testExp;
    uBasicNGSExperiment testExp;
    EXPECT_NO_THROW(testExp.loadWithParser("../data/BED/test.bed","BED" ));
}
TEST(TestExpNOCASE, LOADWIG){

   // uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>, uGenericNGS>  testExp;
     uBasicNGSExperiment testExp;
    EXPECT_NO_THROW(testExp.loadWithParser("../data/WIG/correct.wig","WIG" ));
}
TEST(TestExpNOCASE, LOADSAM){

   // uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>, uGenericNGS>  testExp;
     uBasicNGSExperiment testExp;
    EXPECT_NO_THROW(testExp.loadWithParser("../data/SAM/fiveCountValid.sam","SAM" ));
}

TEST(TestExpNOCASE, LOADBEDCOUNT){

   // uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>, uGenericNGS>  testExp;
     uBasicNGSExperiment testExp;
    ASSERT_NO_THROW(testExp.loadWithParser("../data/BED/test.bed","BED" ));
    EXPECT_EQ(testExp.count(), 4);

    uBasicNGSChrom tempChr = testExp.getChrom("chr3");

    tempChr.sortSites();
    EXPECT_EQ(  tempChr.getSubset(34511,34541).count() ,1 );
}
TEST(TestExpNOCASE, LOADWIGCOUNT){
    //uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>, uGenericNGS>  testExp;
     uBasicNGSExperiment testExp;
    ASSERT_NO_THROW(testExp.loadWithParser("../data/WIG/correct.wig","WIG" ));
    EXPECT_EQ(testExp.count(), 19);
}



