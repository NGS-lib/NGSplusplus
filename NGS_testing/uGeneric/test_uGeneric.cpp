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

#define STARTCASE1 100
#define ENDCASE1 200

/**< Testing our parent unitary uBasicNGS */

/**< Constructor tests */
TEST(uBasicNGSTest, DefaultCTr){
 uBasicNGS uTest;
 EXPECT_EQ(0, uTest.getStart());
 EXPECT_EQ(0, uTest.getEnd());
 EXPECT_EQ("",uTest.getChr());
 EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
}
TEST(uBasicNGSTest, 3CTR){
 uBasicNGS uTest("chr1", STARTCASE1, ENDCASE1);
 EXPECT_EQ(STARTCASE1,uTest.getStart());
 EXPECT_EQ(ENDCASE1, uTest.getEnd());
 EXPECT_EQ("chr1", uTest.getChr());
 EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
}
TEST(uBasicNGSTest, 4CTR){
 uBasicNGS uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::REVERSE);
 EXPECT_EQ(STARTCASE1,uTest.getStart());
 EXPECT_EQ(ENDCASE1, uTest.getEnd());
 EXPECT_EQ("chr1", uTest.getChr());
 EXPECT_EQ(StrandDir::REVERSE,uTest.getStrand());
}
TEST(uBasicNGSTest, 3CTR1SCORE){
 uBasicNGS uTest("chr1", STARTCASE1, ENDCASE1, 0.4f);
 EXPECT_EQ(STARTCASE1,uTest.getStart());
 EXPECT_EQ(ENDCASE1, uTest.getEnd());
 EXPECT_EQ("chr1", uTest.getChr());
 EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
 EXPECT_EQ(0.4f,uTest.getScore());
}

TEST(uBasicNGSTest, 4CTR1SCORE){
 uBasicNGS uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::REVERSE, 0.4);
 EXPECT_EQ(STARTCASE1,uTest.getStart());
 EXPECT_EQ(ENDCASE1, uTest.getEnd());
 EXPECT_EQ("chr1", uTest.getChr());
 EXPECT_EQ(StrandDir::REVERSE,uTest.getStrand());
 EXPECT_EQ(0.4f,uTest.getScore());
}
TEST(uBasicNGSTest, INVCTR){
 EXPECT_ANY_THROW(uBasicNGS uTestInv("chr1",ENDCASE1,STARTCASE1));
}
TEST(uBasicNGSTest, NOCHRCTR){
 EXPECT_NO_THROW(uBasicNGS uTestEmpty("", STARTCASE1, ENDCASE1));
}

/**< End Constructor */

/**< Copy Constructor and implicit conversion */
TEST(uBasicNGSTest, CPYOPP){
 uBasicNGS uTest("chr1", 100, 200,StrandDir::REVERSE, 0.4);

 uBasicNGS uCopy = uTest;
 EXPECT_EQ(STARTCASE1,uCopy.getStart());
 EXPECT_EQ(ENDCASE1, uCopy.getEnd());
 EXPECT_EQ("chr1", uCopy.getChr());
 EXPECT_EQ(StrandDir::REVERSE,uCopy.getStrand());
 EXPECT_EQ(0.4f,uCopy.getScore());

 uTest.setScore(0.5f,2);
 uCopy = uTest;
 EXPECT_EQ(0.5f,uCopy.getScore(2));
}

TEST(uBasicNGSTest, CPYCTR){
 uBasicNGS uTest("chr1", 100, 200,StrandDir::REVERSE, 0.4);
 uTest.setScore(0.5,2);
 uBasicNGS uCopy(uTest);
 EXPECT_EQ(STARTCASE1,uCopy.getStart());
 EXPECT_EQ(ENDCASE1, uCopy.getEnd());
 EXPECT_EQ("chr1", uCopy.getChr());
 EXPECT_EQ(StrandDir::REVERSE,uCopy.getStrand());
 EXPECT_FLOAT_EQ(0.4,uCopy.getScore());
 EXPECT_FLOAT_EQ(0.5,uCopy.getScore(2));
}

TEST(uBasicNGSTest, TAGSCTR){
 uTags uTest("chr1", 100, 200,StrandDir::REVERSE, 0.4);
 uTest.setScore(0.5,2);
 uBasicNGS uCopy(uTest);
 EXPECT_EQ(STARTCASE1,uCopy.getStart());
 EXPECT_EQ(ENDCASE1, uCopy.getEnd());
 EXPECT_EQ("chr1", uCopy.getChr());
 EXPECT_EQ(StrandDir::REVERSE,uCopy.getStrand());
 EXPECT_FLOAT_EQ(0.4,uCopy.getScore());
 EXPECT_FLOAT_EQ(0.5,uCopy.getScore(2));
}

TEST(uBasicNGSTest, REGCTR){
 uRegion uTest("chr1", 100, 200,StrandDir::REVERSE, 0.4);
 uTest.setScore(0.5f,2);
 uBasicNGS uCopy(uTest);
 EXPECT_EQ(STARTCASE1,uCopy.getStart());
 EXPECT_EQ(ENDCASE1, uCopy.getEnd());
 EXPECT_EQ("chr1", uCopy.getChr());
 EXPECT_EQ(StrandDir::REVERSE,uCopy.getStrand());
 EXPECT_FLOAT_EQ(0.4f,uCopy.getScore());
 EXPECT_FLOAT_EQ(0.5f,uCopy.getScore(2));
}
/**< End implicit constructor */


/**< Test Common inherited functions */
TEST(uBasicNGSTest, GETSETSCORE)
    {
        uBasicNGS uTest("chr1", 100, 200,StrandDir::FORWARD, 0.4f);

        EXPECT_ANY_THROW(uTest.getScore(2));
        EXPECT_NO_THROW(uTest.setScore(0.2,1));
        EXPECT_FLOAT_EQ(0.2,uTest.getScore(1));

    }

TEST(uBasicNGSTest, GETSET){
    {
        uBasicNGS uTest("chr1", 100, 200,StrandDir::REVERSE);
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
TEST(uBasicNGSTest, GETSETINV)
{
    uBasicNGS uTest("chr1", 150, 200);

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

    uBasicNGS uTestnew("chr1", 100, 200);
    uBasicNGS uEmpty;

    EXPECT_EQ(101, uTestnew.getLenght());
    /**< Fragment always covers a position. */
    EXPECT_EQ(1, uEmpty.getLenght());
}

TEST(uBasicNGSTest,EXTEND){
    uBasicNGS uTest("chr1", 100, 200);
    uTest.extendSite(100);
    EXPECT_EQ(0, uTest.getStart());
    EXPECT_EQ(300, uTest.getEnd());
    uTest.extendSite(50);
    EXPECT_EQ(0, uTest.getStart());
    EXPECT_EQ(350, uTest.getEnd());
    uBasicNGS uTestExtend("chr1", 500, 600);
    uTestExtend.extendSite(0, 200);
    EXPECT_EQ(500, uTestExtend.getStart());
    EXPECT_EQ(800, uTestExtend.getEnd());
    EXPECT_ANY_THROW(uTest.extendSite(-100));
    EXPECT_ANY_THROW(uTest.extendSite(-100,100));
    EXPECT_ANY_THROW(uTest.extendSite(100,-100));

}
TEST(uBasicNGSTest,TRIM){
    uBasicNGS uTest("chr1", 100, 200);
    uTest.trimSites(50);
    EXPECT_EQ(150, uTest.getStart());
    EXPECT_EQ(150, uTest.getEnd());
    EXPECT_ANY_THROW(uTest.trimSites(50));
    EXPECT_ANY_THROW(uTest.trimSites(-50));

    uBasicNGS uTestTrim2("chr1", 200, 600);
    uTestTrim2.trimSites(50,100);
    EXPECT_EQ(250, uTestTrim2.getStart());
    EXPECT_EQ(500, uTestTrim2.getEnd());
    EXPECT_ANY_THROW(uTest.trimSites(-100,100));
    EXPECT_ANY_THROW(uTest.trimSites(100,-100));
}

TEST(uBasicNGSTest, OVERLAP){

    uBasicNGS uTest("chr1", 100, 200);
    uBasicNGS uTestSame("chr1",100, 200);
    uBasicNGS uTestOverlap1R("chr1", 200, 201);
    uBasicNGS uTestOverlapNot("chr1", 300, 305);
    uBasicNGS uTestOverlapiL("chr1", 99, 100);
    uBasicNGS uTestOverlapDifChr("chr2", 100, 200);
    uTags uTestOverlapPoly("chr1", 100, 200);
    uBasicNGS uTestEmpty;
    uRegion uTestOverlapPolyempty;
    EXPECT_TRUE( uTest.doesOverlap(uTestSame));
    EXPECT_TRUE( uTest.doesOverlap(uTestOverlap1R));
    EXPECT_TRUE( uTest.doesOverlap(uTestOverlapiL));
    EXPECT_TRUE( uTest.doesOverlap(uTestOverlapPoly));
    EXPECT_TRUE(uTestEmpty.doesOverlap(uTestEmpty));

    EXPECT_FALSE (uTest.doesOverlap(uTestEmpty));
    EXPECT_FALSE (uTest.doesOverlap(uTestOverlapPolyempty));
    EXPECT_FALSE (uTest.doesOverlap(uTestOverlapNot));
    EXPECT_FALSE (uTest.doesOverlap(uTestOverlapDifChr));
}

TEST(uBasicNGSTest, DIVIDEINTOBIN){

    uBasicNGS uTest("chr1", 100, 119);

    EXPECT_ANY_THROW(vector<uBasicNGS> ourVector= uTest.divideIntoNBin(3));

    vector<uBasicNGS> TestVector= uTest.divideIntoNBin(3, SplitType::IGNORE);

    EXPECT_EQ( (int)TestVector.size(),3);

    for(uBasicNGS x : TestVector)
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

TEST(uBasicNGSTest, DIVIDEINTOBINOFSIZESTRICT){

    uBasicNGS uTest("chr1", 100, 119);

    EXPECT_ANY_THROW(vector<uBasicNGS> ourVector= uTest.divideIntoBinofSize(7));

    vector<uBasicNGS> ourVector= uTest.divideIntoBinofSize(5);
    for(uBasicNGS x : ourVector)
        EXPECT_EQ( x.getLenght(), 5);
}

TEST(uBasicNGSTest, DIVIDEINTOBINOFSIZEIGNORE){

    uBasicNGS uTest("chr1", 100, 119);
    vector<uBasicNGS> TestVector= uTest.divideIntoBinofSize(7, SplitType::IGNORE);

    EXPECT_EQ( (int)TestVector.size(),2);

    for(uBasicNGS x : TestVector)
        EXPECT_EQ( x.getLenght(), 7);

    EXPECT_EQ( TestVector.at(0).getStart(), 100);
    EXPECT_EQ( TestVector.at(0).getEnd(), 106);
    EXPECT_EQ( TestVector.at(1).getStart(), 107);
    EXPECT_EQ( TestVector.at(1).getEnd(), 113);

}
TEST(uBasicNGSTest, DIVIDEINTOBINOFSIZEADD){
   uBasicNGS uTest("chr1", 100, 119);

    vector<uBasicNGS> TestVector= uTest.divideIntoBinofSize(7, SplitType::ADD);
    EXPECT_EQ( (int)TestVector.size(), 3);
    EXPECT_EQ( TestVector.at(2).getLenght(), 6);
    EXPECT_EQ( TestVector.at(2).getStart(), 114);
    EXPECT_EQ( TestVector.at(2).getEnd(), 119);

    TestVector= uTest.divideIntoBinofSize(10, SplitType::ADD);
    EXPECT_EQ( (int)TestVector.size(), 2);

}
TEST(uBasicNGSTest, DIVIDEINTOBINOFSIZEEXTEND){
   uBasicNGS uTest("chr1", 100, 119);
    vector<uBasicNGS> TestVector= uTest.divideIntoBinofSize(7, SplitType::EXTEND);
    EXPECT_EQ( (int)TestVector.size(), 2);
    EXPECT_EQ( TestVector.at(1).getLenght(), 13);
    EXPECT_EQ( TestVector.at(1).getStart(), 107);
    EXPECT_EQ( TestVector.at(1).getEnd(), 119);
}
