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


