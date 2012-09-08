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

using namespace std;

    #define SOMENUMBER 102343

class uItem : public uGenericNGS{

    int chr;

};

class ourDerivedClass : public uGenericNGSChrom<uGenericNGS> {

public:
  ourDerivedClass():uGenericNGSChrom(){};
  ourDerivedClass(std::string chrom):uGenericNGSChrom(chrom){};
  ourDerivedClass(std::string chrom, long int size):uGenericNGSChrom(chrom,size){};
};

class ourDerivedExperiment : public uGenericNGSExperiment<ourDerivedClass, uGenericNGS> {

};



TEST(ExperimentDivide, DIVIDECHROMINTONBIN){

    ourDerivedExperiment ourExp;
    ourExp.addSite(uGenericNGS("chr21", 100, 199));
    ourExp.addSite(uGenericNGS("chr22", 100, 199));
    ourExp.addSite(uGenericNGS("chr23", 100, 199));
    ourExp.addSite(uGenericNGS("chr24", 100, 199));
    ourExp.divideItemsIntoNBins(4);

    EXPECT_EQ(ourExp.count(),16);
}
TEST(ExperimentDivide, DIVIDECHROMINTOBINOFSIZE){


    ourDerivedExperiment ourExp;

    ourExp.addSite(uGenericNGS("chr21", 100, 199));
    ourExp.addSite(uGenericNGS("chr22", 100, 199));
    ourExp.addSite(uGenericNGS("chr23", 100, 199));
    ourExp.addSite(uGenericNGS("chr24", 100, 199));
    ourExp.divideItemsIntoBinofSize(50);
    EXPECT_EQ(ourExp.count(),8);

}



/**< Testing our base container FUNCTIONS */
/**< Not that we should never implement this without a derived child, so as such the construction is not accesible */
TEST(uGenericNGSChromTest, CTR){

    EXPECT_NO_THROW(ourDerivedClass uChromTest);
    EXPECT_NO_THROW(ourDerivedClass uChromTest2("chr1"));
    EXPECT_NO_THROW(ourDerivedClass uChromTest3(""));
}
TEST(uGenericNGSChromTest, FUNCTIONS){


    ourDerivedClass uChromEmpty;
    ourDerivedClass uChromChr1("chr1");
    ourDerivedClass uChromNull("");
    ourDerivedClass uChromBig("chrBig", SOMENUMBER);

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

    EXPECT_ANY_THROW(factory::makeNGSfromTabString("chr2 400 300"));

    uGenericNGS ourTest=factory::makeNGSfromTabString(("chr2 200 300"));

    EXPECT_EQ(ourTest.getChr(), "chr2");
    EXPECT_EQ(ourTest.getStart(), 200);
    EXPECT_EQ(ourTest.getEnd(), 300);
}

TEST(writeTest, uGenecExpTest){

   // ourDerivedExperiment testExp;
   //uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>,uGenericNGS> testExp;
    uTagsExperiment textExp;
    //testExp.addSite(uGenericNGS("chr21", 300, 800));
    textExp.writeAsBedFile(cout);


}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return  RUN_ALL_TESTS();;
}
