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

using namespace std;

    #define SOMENUMBER 102343

using namespace NGS;


TEST(ExperimentDivide, DIVIDECHROMINTONBIN){

    uBasicNGSExperiment ourExp;
    ourExp.addSite(uBasicNGS("chr21", 100, 199));
    ourExp.addSite(uBasicNGS("chr22", 100, 199));
    ourExp.addSite(uBasicNGS("chr23", 100, 199));
    ourExp.addSite(uBasicNGS("chr24", 100, 199));
    ourExp.divideItemsIntoNBins(4);

    EXPECT_EQ(ourExp.count(),16);
}
TEST(ExperimentDivide, DIVIDECHROMINTOBINOFSIZE){


    uBasicNGSExperiment ourExp;

    ourExp.addSite(uBasicNGS("chr21", 100, 199));
    ourExp.addSite(uBasicNGS("chr22", 100, 199));
    ourExp.addSite(uBasicNGS("chr23", 100, 199));
    ourExp.addSite(uBasicNGS("chr24", 100, 199));
    ourExp.divideItemsIntoBinofSize(50);
    EXPECT_EQ(ourExp.count(),8);

}



/**< Testing our base container FUNCTIONS */
/**< Not that we should never implement this without a derived child, so as such the construction is not accesible */
TEST(uBasicNGSChromTest, CTR){

    EXPECT_NO_THROW(uBasicNGSChrom uChromTest);
    EXPECT_NO_THROW(uBasicNGSChrom uChromTest2("chr1"));
    EXPECT_NO_THROW(uBasicNGSChrom uChromTest3(""));
}
TEST(uBasicNGSChromTest, FUNCTIONS){


    uBasicNGSChrom uChromEmpty;
    uBasicNGSChrom uChromChr1("chr1");
    uBasicNGSChrom uChromNull("");
    uBasicNGSChrom uChromBig("chrBig", SOMENUMBER);

    EXPECT_EQ(uChromEmpty.getChr() , uChromNull.getChr());

    EXPECT_EQ("chr1", uChromChr1.getChr());
    EXPECT_EQ(0, uChromChr1.count());
    EXPECT_EQ(0, uChromChr1.getChromSize());
    EXPECT_ANY_THROW(uChromChr1.setChromSize(-1));
    EXPECT_NO_THROW(uChromChr1.setChromSize(SOMENUMBER));
    EXPECT_EQ(uChromBig.getChromSize(),uChromChr1.getChromSize());

    uChromBig.setChr("chr3");
    EXPECT_EQ("chr3", uChromBig.getChr());

    EXPECT_ANY_THROW(uChromBig.getSite(0));
}


TEST(factoryTest, uGenericTest){
   // std::string myline("chr2 400 300");
  //  EXPECT_ANY_THROW(factory::makeNGSfromTabString(myline));

    uBasicNGS ourTest=factory::makeNGSfromTabString<uBasicNGS>(("chr2 200 300"));

    EXPECT_EQ(ourTest.getChr(), "chr2");
    EXPECT_EQ(ourTest.getStart(), 200);
    EXPECT_EQ(ourTest.getEnd(), 300);
}

TEST(writeTest, uGenecExpTest){

   // uBasicNGSExperiment testExp;
   //uBasicNGSExperiment<uBasicNGSChrom<uBasicNGS>,uBasicNGS> testExp;
    uTagsExperiment textExp;
    //testExp.addSite(uBasicNGS("chr21", 300, 800));
    textExp.writeAsBedFile(cout);


}

using namespace testing;
int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
 return  RUN_ALL_TESTS();;
}
