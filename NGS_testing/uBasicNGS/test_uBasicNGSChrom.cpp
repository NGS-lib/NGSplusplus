
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
// To use a test fixture, derive a class from testing::Test.à
#define CHROMDIVIDESIZE 500
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

class ChromSubset : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {
    uChromSubsetTest.setChr("chr1");
    uChromSubsetTest.addData(uBasicNGS("chr1", 300, 500));
    uChromSubsetTest.addData(uBasicNGS("chr1", 100, 199));
    uChromSubsetTest.addData(uBasicNGS("chr1", 100, 299));
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
    uChromTestOverlap.addData(uBasicNGS("chr1", 200, 800));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 800));
    uChromTestOverlap.addData(uBasicNGS("chr1", 300, 800));
    EXPECT_FALSE(uChromTestOverlap.isSorted());
}

