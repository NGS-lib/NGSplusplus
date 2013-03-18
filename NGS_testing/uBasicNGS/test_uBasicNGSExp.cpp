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

class TestExp : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {
    uExpTest.addData(uBasicNGS("chr1", 300, 500, 0.1f));
    uExpTest.addData(uBasicNGS("chr1", 100, 199,0.5f));
    uExpTest.addData(uBasicNGS("chr1", 100, 299,0.3f));
  }
uBasicNGSExperiment uExpTest;
};

TEST_F(TestExp, SORTEXP){
    /**< Unsorted, fail */
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites();
    EXPECT_TRUE(uExpTest.isSorted());
    uExpTest.addData(uBasicNGS("chr1", 2, 50));
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites();
    EXPECT_TRUE(uExpTest.isSorted());
}
TEST_F(TestExp, SORTCUSTOM){

    auto comp=[](const uBasicNGS & item1, const uBasicNGS & item2){
             return item1.getScore() < item2.getScore();
    };
    /**< Unsorted, fail */
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites(comp);
    EXPECT_TRUE(uExpTest.isSorted());
    uExpTest.addData(uBasicNGS("chr1", 2, 50, 0.05f));
    uExpTest.addData(uBasicNGS("chr1", 2, 50, 0.7f));
    EXPECT_FALSE(uExpTest.isSorted());
    uExpTest.sortSites();
    EXPECT_TRUE(uExpTest.isSorted());
}


TEST_F(TestExp, GETSUBSET){

    uExpTest.addData(uBasicNGS("chr1", 2, 50, 0.05f));
    uExpTest.addData(uBasicNGS("chr1", 5, 25, 0.7f));
    uExpTest.sortSites();
    auto resultChrom=uExpTest.getSubset("chr1",99,200);
    EXPECT_EQ(resultChrom.count(), 2);

}


TEST_F(TestExp, GETSUBSETLARGE){

    uExpTest.addData(uBasicNGS("chr1", 180000050, 180500000 ));
    uExpTest.addData(uBasicNGS("chr1", 180050000, 180500000));
    uExpTest.addData(uBasicNGS("chr1", 180550000, 180900000));
    uExpTest.addData(uBasicNGS("chr1", 5, 25, 0.7f));
    uExpTest.addData(uBasicNGS("chr1", 5, 25, 0.7f));

    uExpTest.sortSites();
    auto resultChrom=uExpTest.getSubset("chr1",180000000,181000000);
    EXPECT_EQ(resultChrom.count(), 3);

}

TEST_F(TestExp, GETSUBSETCUSTOM){

  auto comp=[](const uBasicNGS & item1, const uBasicNGS & item2){
             return item1.getScore() < item2.getScore();
    };

    uExpTest.addData(uBasicNGS("chr1", 2, 50, 0.05f));
    uExpTest.addData(uBasicNGS("chr1", 5, 25, 0.7f));
    auto getScoreFunct=std::mem_fn((float(uBasicNGS::*)()const)&uBasicNGS::getScore);

    uExpTest.sortSites(comp,getScoreFunct,getScoreFunct);
    auto resultChrom=uExpTest.getSubset("chr1",0.4f,0.9f);
    EXPECT_EQ(resultChrom.count(), 3);
}


TEST(TestExpNOCASE, LOADBED){

   // uBasicNGSExperiment<uBasicNGSChrom<uBasicNGS>, uBasicNGS>  testExp;
    uBasicNGSExperiment testExp;
    EXPECT_NO_THROW(testExp.loadWithParser("../data/BED/test.bed","BED" ));
}
TEST(TestExpNOCASE, LOADWIG){

   // uBasicNGSExperiment<uBasicNGSChrom<uBasicNGS>, uBasicNGS>  testExp;
     uBasicNGSExperiment testExp;
    EXPECT_NO_THROW(testExp.loadWithParser("../data/WIG/correct.wig","WIG" ));
}
TEST(TestExpNOCASE, LOADSAM){

   // uBasicNGSExperiment<uBasicNGSChrom<uBasicNGS>, uBasicNGS>  testExp;
     uBasicNGSExperiment testExp;
    EXPECT_NO_THROW(testExp.loadWithParser("../data/SAM/fiveCountValid.sam","SAM" ));
}

TEST(TestExpNOCASE, LOADBEDCOUNT){

   // uBasicNGSExperiment<uBasicNGSChrom<uBasicNGS>, uBasicNGS>  testExp;
    uBasicNGSExperiment testExp;
    ASSERT_NO_THROW(testExp.loadWithParser("../data/BED/test.bed","BED" ));
    EXPECT_EQ(testExp.count(), 4);

    uBasicNGSChrom tempChr = testExp.getChrom("chr3");

    tempChr.sortSites();
    EXPECT_EQ(  tempChr.getSubset(34511,34541).count() ,1 );
}
TEST(TestExpNOCASE, LOADWIGCOUNT){
    //uBasicNGSExperiment<uBasicNGSChrom<uBasicNGS>, uBasicNGS>  testExp;
     uBasicNGSExperiment testExp;
    ASSERT_NO_THROW(testExp.loadWithParser("../data/WIG/correct.wig","WIG" ));
    EXPECT_EQ(testExp.count(), 19);
}


TEST(uBasicNGSsExpTest_copyCtr, NORMAL){
    ASSERT_TRUE(false);
}




