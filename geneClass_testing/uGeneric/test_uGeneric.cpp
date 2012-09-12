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
using namespace NGS;
/**< Testing our parent unitary uGenericNGS */
TEST(uGenericNGSTest, DefaultCTr){
 uGenericNGS uTest;
 EXPECT_EQ(0, uTest.getStart());
 EXPECT_EQ(0, uTest.getEnd());
 EXPECT_EQ("",uTest.getChr());
 EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
}

TEST(uGenericNGSTest, 3CTR){
 uGenericNGS uTest("chr1", 100, 200);
 EXPECT_EQ(100,uTest.getStart());
 EXPECT_EQ(200, uTest.getEnd());
 EXPECT_EQ("chr1", uTest.getChr());
}
TEST(uGenericNGSTest, INVCTR){
 EXPECT_ANY_THROW(uGenericNGS uTestInv("chr1",200,100));
}
TEST(uGenericNGSTest, NOCHRCTR){
 EXPECT_NO_THROW(uGenericNGS uTestEmpty("", 100, 200));
}
TEST(uGenericNGSTest, GETSET){

    {
        uGenericNGS uTest("chr1", 100, 200,StrandDir::REVERSE);
        /**< Illegal Start */
        EXPECT_ANY_THROW(uTest.setStart(-10));
        EXPECT_ANY_THROW(uTest.setStart(250));
        EXPECT_EQ(100,uTest.getStart());
        uTest.setStart(150);
        EXPECT_EQ(150, uTest.getStart());
        uTest.setStart(200);
        EXPECT_EQ(200, uTest.getStart());
        uTest.setStart(150);

        EXPECT_EQ(StrandDir::REVERSE,uTest.getStrand());
        uTest.setStrand(StrandDir::FORWARD);
        EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());

        /**< Strand char set */
        uTest.setStrand('-');
        EXPECT_EQ(StrandDir::REVERSE,uTest.getStrand());
        EXPECT_THROW(uTest.setStrand('a'), param_throw);
        EXPECT_TRUE(uTest.isReverse());
    }
 }
TEST(uGenericNGSTest, GETSETINV)
{
    uGenericNGS uTest("chr1", 150, 200);

    /**< Illegal end */
    EXPECT_ANY_THROW(uTest.setEnd(-10));
    EXPECT_ANY_THROW(uTest.setEnd(50));
    EXPECT_EQ(200,uTest.getEnd());

    EXPECT_NO_THROW(uTest.setEnd(150));

    EXPECT_EQ(150,uTest.getEnd());
    EXPECT_NO_THROW(uTest.setEnd(225));
    EXPECT_EQ(225,uTest.getEnd());

    /**< Chrom setss */
    EXPECT_EQ("chr1",uTest.getChr());
    uTest.setChr("chr2");
    EXPECT_EQ("chr2",uTest.getChr());

    uGenericNGS uTestnew("chr1", 100, 200);
    uGenericNGS uEmpty;

    EXPECT_EQ(101, uTestnew.getLenght());

    /**< Fragment always covers a position. */
    EXPECT_EQ(1, uEmpty.getLenght());

}

TEST(uGenericNGSTest,EXTEND){
    uGenericNGS uTest("chr1", 100, 200);
    uTest.extendSite(100);
    EXPECT_EQ(0, uTest.getStart());
    EXPECT_EQ(300, uTest.getEnd());
    uTest.extendSite(50);
    EXPECT_EQ(0, uTest.getStart());
    EXPECT_EQ(350, uTest.getEnd());
    uGenericNGS uTestExtend("chr1", 500, 600);
    uTestExtend.extendSite(0, 200);
    EXPECT_EQ(500, uTestExtend.getStart());
    EXPECT_EQ(800, uTestExtend.getEnd());
    EXPECT_ANY_THROW(uTest.extendSite(-100));
    EXPECT_ANY_THROW(uTest.extendSite(-100,100));
    EXPECT_ANY_THROW(uTest.extendSite(100,-100));

}
TEST(uGenericNGSTest,TRIM){
    uGenericNGS uTest("chr1", 100, 200);
    uTest.trimSites(50);
    EXPECT_EQ(150, uTest.getStart());
    EXPECT_EQ(150, uTest.getEnd());
    EXPECT_ANY_THROW(uTest.trimSites(50));
    EXPECT_ANY_THROW(uTest.trimSites(-50));

    uGenericNGS uTestTrim2("chr1", 200, 600);
    uTestTrim2.trimSites(50,100);
    EXPECT_EQ(250, uTestTrim2.getStart());
    EXPECT_EQ(500, uTestTrim2.getEnd());
    EXPECT_ANY_THROW(uTest.trimSites(-100,100));
    EXPECT_ANY_THROW(uTest.trimSites(100,-100));
}

TEST(uGenericNGSTest, OVERLAP){

     uGenericNGS uTest("chr1", 100, 200);
     uGenericNGS uTestSame("chr1",100, 200);
     uGenericNGS uTestOverlap1R("chr1", 200, 201);
     uGenericNGS uTestOverlapNot("chr1", 300, 305);
     uGenericNGS uTestOverlapiL("chr1", 99, 100);
     uGenericNGS uTestOverlapDifChr("chr2", 100, 200);
     uTags uTestOverlapPoly("chr1", 100, 200);
     uGenericNGS uTestEmpty;
     uRegion uTestOverlapPolyempty;
    std::cerr << "Overlap 1" <<std::endl;
    EXPECT_TRUE( uTest.doesOverlap(uTestSame));
    EXPECT_TRUE( uTest.doesOverlap(uTestOverlap1R));
    EXPECT_TRUE( uTest.doesOverlap(uTestOverlapiL));
    EXPECT_TRUE( uTest.doesOverlap(uTestOverlapPoly));
    EXPECT_TRUE(uTestEmpty.doesOverlap(uTestEmpty));

    std::cerr << "Overlap 2" <<std::endl;
    EXPECT_FALSE (uTest.doesOverlap(uTestEmpty));
    EXPECT_FALSE (uTest.doesOverlap(uTestOverlapPolyempty));
    EXPECT_FALSE (uTest.doesOverlap(uTestOverlapNot));
    EXPECT_FALSE (uTest.doesOverlap(uTestOverlapDifChr));
}

TEST(uGenericNGSTest, DIVIDEINTOBIN){

    uGenericNGS uTest("chr1", 100, 119);

    EXPECT_ANY_THROW(vector<uGenericNGS> ourVector= uTest.divideIntoNBin(3));

    vector<uGenericNGS> TestVector= uTest.divideIntoNBin(3, SplitType::IGNORE);

    EXPECT_EQ( (int)TestVector.size(),3);

    for(uGenericNGS x : TestVector)
        EXPECT_EQ( x.getLenght(), 6);

    EXPECT_EQ( TestVector.at(0).getStart(), 100);
    EXPECT_EQ( TestVector.at(0).getEnd(), 105);
    EXPECT_EQ( TestVector.at(1).getStart(), 106);
    EXPECT_EQ( TestVector.at(1).getEnd(), 111);
    EXPECT_EQ( TestVector.at(2).getStart(), 112);
    EXPECT_EQ( TestVector.at(2).getEnd(), 117);

    TestVector= uTest.divideIntoNBin(3, SplitType::ADD);
    EXPECT_EQ( (int)TestVector.size(), 4);
    EXPECT_EQ( TestVector.at(3).getLenght(), 2);
    EXPECT_EQ( TestVector.at(3).getStart(), 118);
    EXPECT_EQ( TestVector.at(3).getEnd(), 119);

    TestVector= uTest.divideIntoNBin(4, SplitType::ADD);
    EXPECT_EQ( (int)TestVector.size(), 4);

    TestVector= uTest.divideIntoNBin(3, SplitType::EXTEND);
    EXPECT_EQ( (int)TestVector.size(), 3);
    EXPECT_EQ( TestVector.at(2).getLenght(), 8);
    EXPECT_EQ( TestVector.at(2).getStart(), 112);
    EXPECT_EQ( TestVector.at(2).getEnd(), 119);
}

TEST(uGenericNGSTest, DIVIDEINTOBINOFSIZESTRICT){

    uGenericNGS uTest("chr1", 100, 119);

    EXPECT_ANY_THROW(vector<uGenericNGS> ourVector= uTest.divideIntoBinofSize(7));

    vector<uGenericNGS> ourVector= uTest.divideIntoBinofSize(5);
    for(uGenericNGS x : ourVector)
        EXPECT_EQ( x.getLenght(), 5);
}

TEST(uGenericNGSTest, DIVIDEINTOBINOFSIZEIGNORE){

    uGenericNGS uTest("chr1", 100, 119);
    vector<uGenericNGS> TestVector= uTest.divideIntoBinofSize(7, SplitType::IGNORE);

    EXPECT_EQ( (int)TestVector.size(),2);

    for(uGenericNGS x : TestVector)
        EXPECT_EQ( x.getLenght(), 7);

    EXPECT_EQ( TestVector.at(0).getStart(), 100);
    EXPECT_EQ( TestVector.at(0).getEnd(), 106);
    EXPECT_EQ( TestVector.at(1).getStart(), 107);
    EXPECT_EQ( TestVector.at(1).getEnd(), 113);

}
TEST(uGenericNGSTest, DIVIDEINTOBINOFSIZEADD){
   uGenericNGS uTest("chr1", 100, 119);

    vector<uGenericNGS> TestVector= uTest.divideIntoBinofSize(7, SplitType::ADD);
    EXPECT_EQ( (int)TestVector.size(), 3);
    EXPECT_EQ( TestVector.at(2).getLenght(), 6);
    EXPECT_EQ( TestVector.at(2).getStart(), 114);
    EXPECT_EQ( TestVector.at(2).getEnd(), 119);

    TestVector= uTest.divideIntoBinofSize(10, SplitType::ADD);
    EXPECT_EQ( (int)TestVector.size(), 2);

}
TEST(uGenericNGSTest, DIVIDEINTOBINOFSIZEEXTEND){
   uGenericNGS uTest("chr1", 100, 119);

    vector<uGenericNGS> TestVector= uTest.divideIntoBinofSize(7, SplitType::EXTEND);
    EXPECT_EQ( (int)TestVector.size(), 2);
    EXPECT_EQ( TestVector.at(1).getLenght(), 13);
    EXPECT_EQ( TestVector.at(1).getStart(), 107);
    EXPECT_EQ( TestVector.at(1).getEnd(), 119);

}
