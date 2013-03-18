
/**< Test Common inherited functions */
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
#include <gtest/gtest.h>

using namespace std;
using namespace NGS;

#define SOMENUMBER 102343
#define STARTCASE1 100
#define ENDCASE1 200










/**<


Test Common inherited functions


NEVER MODIFY THESE TESTS DIRECTLY. They are
a carbon copy of the uBasic Tests and are included to validate no performance
change in inherited versions.


*/




/*
 * Setters/Getters testing - start/end (positions)
 *	Valid cases:
 *		SETGETSTART
 *		SETGETEND
 *	Invalid cases:
 *		SETSTARTILLEGAL
 *		SETSENDILLEGAL
 */

TEST(uRegionNGSTestHerit, SETGETSTART){
        uRegion uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::REVERSE);
        EXPECT_EQ(STARTCASE1,uTest.getStart());
        uTest.setStart(150);
        EXPECT_EQ(150, uTest.getStart());
        uTest.setStart(ENDCASE1);
        EXPECT_EQ(ENDCASE1, uTest.getStart());
 }

TEST(uRegionNGSTestHerit, SETGETEND){
        uRegion uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::REVERSE);
        EXPECT_EQ(ENDCASE1,uTest.getEnd());
        uTest.setEnd(250);
        EXPECT_EQ(250, uTest.getEnd());
        uTest.setEnd(100);
        EXPECT_EQ(100, uTest.getEnd());
 }

TEST(uRegionNGSTestHerit, SETSTARTILLEGAL){
        uRegion uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::REVERSE);
        /**< Illegal Start */
        EXPECT_ANY_THROW(uTest.setStart(-10));
        EXPECT_ANY_THROW(uTest.setStart(250));
 }

 TEST(uRegionNGSTestHerit, SETSENDILLEGAL){
        uRegion uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::REVERSE);
        /**< Illegal Start */
        EXPECT_ANY_THROW(uTest.setEnd(50));
        EXPECT_ANY_THROW(uTest.setEnd(-20));
 }

/**< TODO: Add test for StartEnd() */

/*
 * Setters/Getters testing - strand
 *	Valid cases:
 *		SETGETSTRAND
 *		SETGETSTRANDCHAR
 *	Invalid case:
 *		SETSTRANDFAIL
 */

TEST(uRegionNGSTestHerit, SETGETSTRAND){
        uRegion uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::FORWARD);
        EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
        uTest.setStrand(StrandDir::REVERSE);
        EXPECT_EQ(StrandDir::REVERSE,uTest.getStrand());
        EXPECT_TRUE(uTest.isReverse());
 }

TEST(uRegionNGSTestHerit, SETGETSTRANDCHAR){
        uRegion uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::FORWARD);
        uTest.setStrand('+');
        EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
        uTest.setStrand('-');
        EXPECT_EQ(StrandDir::REVERSE,uTest.getStrand());
}

TEST(uRegionNGSTestHerit, SETSTRANDFAIL){
       uRegion uTest("chr1", STARTCASE1, ENDCASE1,StrandDir::FORWARD);
       EXPECT_THROW(uTest.setStrand('a'), param_throw);
}

/*
 * Setters/Getters testing - chr
 *	Valid case:
 *		GETSETCHR
 *	Invalid case:
 */

TEST(uRegionNGSTestHerit, GETSETCHR)
{
    uRegion uTest("chr1", 150, 200);
    EXPECT_EQ("chr1",uTest.getChr());
    uTest.setChr("chr2");
    EXPECT_EQ("chr2",uTest.getChr());
    uRegion uEmpty;
    EXPECT_EQ("",uEmpty.getChr());
}

/*
 * Setters/Getters testing - length
 *	Valid case:
 *		GETLENGHT
 *	Invalid case:
 */

TEST(uRegionNGSTestHerit, GETLENGHT)
{
    uRegion uTest("chr1", 150, 200);
    EXPECT_EQ(51, uTest.getLenght());
    /**< Fragment always covers a position. */
    uRegion uEmpty;
    EXPECT_EQ(1, uEmpty.getLenght());
    uTest.setStart(200);
    EXPECT_EQ(1, uTest.getLenght());
}

/*
 * Setters/Getters testing - score
 *	Valid cases:
 *		GETSETSCORE
 *		GETSETSCOREVECTOR
 *	Invalid case:
 */

TEST(uRegionNGSTestHerit, GETSETSCORE)
{
    uRegion uTest("chr1", 100, 200,StrandDir::FORWARD, 0.4f);
    EXPECT_FLOAT_EQ(0.4f,uTest.getScore(0));
    EXPECT_ANY_THROW(uTest.getScore(2));
    EXPECT_NO_THROW(uTest.setScore(0.2,1));
    EXPECT_FLOAT_EQ(0.2,uTest.getScore(1));
}

TEST(uRegionNGSTestHerit, GETSETSCOREVECTOR)
{
   uRegion uTest("chr1", 100, 200,StrandDir::FORWARD);
   uTest.setScoreVector({4.2f,0.4f,-1.03f});
   EXPECT_FLOAT_EQ(4.2f, uTest.getScore());
   EXPECT_EQ(3, (int)(uTest.getScoreVector().size()));
   auto vec=uTest.getScoreVector();
   EXPECT_FLOAT_EQ(4.2f,vec.at(0));
   EXPECT_FLOAT_EQ(0.4f,vec.at(1));
   EXPECT_FLOAT_EQ(-1.03f,vec.at(2));
}

/**< TODO: Is it possible to have an invalid score? */

/*
 * Test for extending functions
 *		void extendSite(int extend);
 *		void extendSite(int extendLeft, int extendRight);
 *
 *	Valid cases:
 *		EXTENDSIMPLE
 *		EXTENDDOUBLE
 *		EXTENDNOTHING
 *	Invalid case:
 *		EXTENDFAIL
 */

TEST(uRegionNGSTestHerit,EXTENDSIMPLE){
    uRegion uTest("chr1", 100, 200);
    uTest.extendSite(100);
    EXPECT_EQ(0, uTest.getStart());
    EXPECT_EQ(300, uTest.getEnd());
    uTest.extendSite(50);
    EXPECT_EQ(0, uTest.getStart());
    EXPECT_EQ(350, uTest.getEnd());
}

TEST(uRegionNGSTestHerit,EXTENDDOUBLE){
    uRegion uTestExtend("chr1", 500, 600);
    uTestExtend.extendSite(0, 200);
    EXPECT_EQ(500, uTestExtend.getStart());
    EXPECT_EQ(800, uTestExtend.getEnd());

}

TEST(uRegionNGSTestHerit,EXTENDNOTHING){
    uRegion uTest("chr1", 500, 600);
    EXPECT_NO_THROW(uTest.extendSite(0));
    EXPECT_EQ(500, uTest.getStart());
    EXPECT_EQ(600, uTest.getEnd());
}

TEST(uRegionNGSTestHerit,EXTENDFAIL){
    uRegion uTest("chr1", 500, 600);
    EXPECT_ANY_THROW(uTest.extendSite(-100));
    EXPECT_ANY_THROW(uTest.extendSite(-100,100));
    EXPECT_ANY_THROW(uTest.extendSite(100,-100));
}

/*
 * Test for trimming functions
 *		void trimSite(int trim);
 *		void trimSite(int trimLeft,int trimRight);
 *
 *	Valid cases:
 *		TRIMSIMPLE
 *		TRIMDOUBLE
 *		TRIMNOTHING
 *	Invalid case:
 *		TRIMILLEGAL
 */

//TODO, verify behavior when trimming beyond 0.
TEST(uRegionNGSTestHerit,TRIMSIMPLE){
    uRegion uTest("chr1", 100, 200);
    uTest.trimSite(50);
    EXPECT_EQ(150, uTest.getStart());
    EXPECT_EQ(150, uTest.getEnd());
    EXPECT_ANY_THROW(uTest.trimSite(50));
    EXPECT_ANY_THROW(uTest.trimSite(-50));
}

TEST(uRegionNGSTestHerit,TRIMDOUBLE){
    uRegion uTest("chr1", 200, 600);
    uTest.trimSite(50,100);
    EXPECT_EQ(250, uTest.getStart());
    EXPECT_EQ(500, uTest.getEnd());
}

TEST(uRegionNGSTestHerit,TRIMILLEGAL){
    uRegion uTestTrim2("chr1", 200, 600);
    EXPECT_ANY_THROW(uTestTrim2.trimSite(-100,100));
    EXPECT_ANY_THROW(uTestTrim2.trimSite(100,-100));
}

TEST(uRegionNGSTestHerit,TRIMNOTHING){
    uRegion uTest("chr1", 500, 600);
    EXPECT_NO_THROW(uTest.trimSite(0));
    EXPECT_EQ(500, uTest.getStart());
    EXPECT_EQ(600, uTest.getEnd());
}

/*
 * Test for overlap function
 *		bool doesOverlap(_SELF_ other,OverlapType type=OverlapType::OVERLAP_PARTIAL) const;
 *
 *	Valid cases:
 *		OVERLAP
 *		OVERLAPDIFFCHR
 *		OVERLAPNOT
 *	Invalid case:
 */

TEST(uRegionNGSTestHerit, OVERLAP){
    uRegion uTest("chr1", 100, 200);
    uRegion uTestSame("chr1",100, 200);
    uRegion uTestOverlap1R("chr1", 200, 201);
    uRegion uTestOverlapiL("chr1", 99, 100);

    EXPECT_TRUE( uTest.doesOverlap(uTestSame));
    EXPECT_TRUE( uTest.doesOverlap(uTestOverlap1R));
    EXPECT_TRUE( uTest.doesOverlap(uTestOverlapiL));
    EXPECT_TRUE( uTest.doesOverlap(uTest));
}

TEST(uRegionNGSTestHerit, OVERLAPDIFFCHR){
    uRegion uTest("chr1", 100, 200);
    uRegion uTestOverlapDifChr("chr2", 100, 200);
    uRegion uTestEmpty;
    EXPECT_FALSE (uTest.doesOverlap(uTestOverlapDifChr));
    EXPECT_FALSE (uTest.doesOverlap(uTestEmpty));
}

TEST(uRegionNGSTestHerit, OVERLAPNOT){
    uRegion uTest("chr1", 100, 200);
    uRegion uTestOverlapNot("chr1", 300, 305);
    uRegion uTestNot2("chr1", 201, 300);
    EXPECT_FALSE (uTest.doesOverlap(uTestOverlapNot));
    EXPECT_FALSE (uTest.doesOverlap(uTestNot2));
}

/*
 * Test for divide into bin functions
 *		bool doesOverlap(_SELF_ other,OverlapType type=OverlapType::OVERLAP_PARTIAL) const;
 *
 *	Valid cases:
 *		DIVIDEINTOBINADD
 *		DIVIDEINTOBINEXTEND
 *		DIVIDEINTOBINIGNORE
 *		DIVIDEINTOBINOFSIZESTRICT
 *		DIVIDEINTOBINOFSIZEIGNORE
 *		DIVIDEINTOBINOFSIZEADD
 *		DIVIDEINTOBINOFSIZEEXTEND
 *	Invalid case:
 */

/**< Need to test for every overload */
TEST(uRegionNGSTestHerit, DIVIDEINTOBINADD){

    uRegion uTest("chr1", 100, 119);

    EXPECT_ANY_THROW(vector<uRegion> ourVector= uTest.divideIntoNBin(3));

    auto TestVector= uTest.divideIntoNBin(3, SplitType::ADD);
    EXPECT_EQ( TestVector.at(3).getLenght(), 2);
    EXPECT_EQ( TestVector.at(3).getStart(), 118);
    EXPECT_EQ( TestVector.at(3).getEnd(), 119);

    EXPECT_EQ( 4,(int)TestVector.size());

}

TEST(uRegionNGSTestHerit, DIVIDEINTOBINEXTEND){
     uRegion uTest("chr1", 100, 119);
     auto TestVector= uTest.divideIntoNBin(3, SplitType::EXTEND);
     EXPECT_EQ( (int)TestVector.size(), 3);
     EXPECT_EQ( TestVector.at(2).getLenght(), 8);
     EXPECT_EQ( TestVector.at(2).getStart(), 112);
     EXPECT_EQ( TestVector.at(2).getEnd(), 119);
}

TEST(uRegionNGSTestHerit, DIVIDEINTOBINIGNORE){

    uRegion uTest("chr1", 100, 119);

    EXPECT_ANY_THROW(vector<uRegion> ourVector= uTest.divideIntoNBin(3));

    vector<uRegion> TestVector= uTest.divideIntoNBin(3, SplitType::IGNORE);

    EXPECT_EQ( (int)TestVector.size(),3);

    /**< Test equal size */
    for(uRegion x : TestVector)
        EXPECT_EQ( x.getLenght(), 6);
    EXPECT_EQ( TestVector.at(0).getStart(), 100);
    EXPECT_EQ( TestVector.at(0).getEnd(), 105);
    EXPECT_EQ( TestVector.at(1).getStart(), 106);
    EXPECT_EQ( TestVector.at(1).getEnd(), 111);
    EXPECT_EQ( TestVector.at(2).getStart(), 112);
    EXPECT_EQ( TestVector.at(2).getEnd(), 117);

}

TEST(uRegionNGSTestHerit, DIVIDEINTOBINOFSIZESTRICT){

    uRegion uTest("chr1", 100, 119);

    EXPECT_ANY_THROW(vector<uRegion> ourVector= uTest.divideIntoBinofSize(7));

    vector<uRegion> ourVector= uTest.divideIntoBinofSize(5);
    for(uRegion x : ourVector)
        EXPECT_EQ( x.getLenght(), 5);
}

//TODO
/**< add extra testing to element, test first, last and middle */
TEST(uRegionNGSTestHerit, DIVIDEINTOBINOFSIZEIGNORE){

    uRegion uTest("chr1", 100, 119);
    vector<uRegion> TestVector= uTest.divideIntoBinofSize(7, SplitType::IGNORE);

    EXPECT_EQ( (int)TestVector.size(),2);

    for(uRegion x : TestVector)
        EXPECT_EQ( x.getLenght(), 7);

    EXPECT_EQ( TestVector.at(0).getStart(), 100);
    EXPECT_EQ( TestVector.at(0).getEnd(), 106);
    EXPECT_EQ( TestVector.at(1).getStart(), 107);
    EXPECT_EQ( TestVector.at(1).getEnd(), 113);
}

TEST(uRegionNGSTestHerit, DIVIDEINTOBINOFSIZEADD){
   uRegion uTest("chr1", 100, 119);

    vector<uRegion> TestVector= uTest.divideIntoBinofSize(7, SplitType::ADD);
    EXPECT_EQ( (int)TestVector.size(), 3);
    EXPECT_EQ( TestVector.at(2).getLenght(), 6);
    EXPECT_EQ( TestVector.at(2).getStart(), 114);
    EXPECT_EQ( TestVector.at(2).getEnd(), 119);

    TestVector= uTest.divideIntoBinofSize(10, SplitType::ADD);
    EXPECT_EQ( (int)TestVector.size(), 2);

}
TEST(uRegionNGSTestHerit, DIVIDEINTOBINOFSIZEEXTEND){
   uRegion uTest("chr1", 100, 119);
    vector<uRegion> TestVector= uTest.divideIntoBinofSize(7, SplitType::EXTEND);
    EXPECT_EQ( (int)TestVector.size(), 2);
    EXPECT_EQ( TestVector.at(1).getLenght(), 13);
    EXPECT_EQ( TestVector.at(1).getStart(), 107);
    EXPECT_EQ( TestVector.at(1).getEnd(), 119);
}

/*
 * Test for divide into bin functions
 *		virtual uToken createToken()const;
 *
 *	Valid case:
 *		GETTOKENTEST
 *	Invalid case:
 */

TEST(uRegionNGSTestHerit, GETTOKENTEST){
    uRegion uTest("chr1", 100, 119);

    uToken ourToken= uTest.createToken();
    EXPECT_EQ("chr1", ourToken.getParam(token_param::CHR));
    EXPECT_EQ("100", ourToken.getParam(token_param::START_POS));
    EXPECT_EQ("119", ourToken.getParam(token_param::END_POS));
    /**< Random Testin of certain Params */
    EXPECT_FALSE(ourToken.isParamSet(token_param::CIGAR));
    EXPECT_FALSE(ourToken.isParamSet(token_param::DENSITY));
    EXPECT_FALSE(ourToken.isParamSet(token_param::SCORE));
}

/**< isEqual */

TEST(uRegionNGSTestHerit, ISEQUAL){
    uRegion uTest("chr1", 100, 119);
    uRegion uTestCopy("chr1", 100, 119);
    uRegion uTestChrChange("chr2", 100, 119);
    uRegion uTestStartChange("chr1", 110, 119);
    uRegion uTestEndChange("chr1", 100, 129);
    uRegion uTestScoreChange("chr1", 100, 119,0.4f);

    EXPECT_TRUE(uTest.isEqual(uTestCopy));
    EXPECT_FALSE(uTest.isEqual(uTestChrChange));
    EXPECT_FALSE(uTest.isEqual(uTestStartChange));
    EXPECT_FALSE(uTest.isEqual(uTestEndChange));
    EXPECT_FALSE(uTest.isEqual(uTestScoreChange));
}

/**< GetCopy */

TEST(uRegionNGSTestHerit, GETCOPY){
    uRegion uTest("chr1", 100, 119);

   EXPECT_TRUE( uTest.isEqual(uTest.getCopy()));
}


/**< Print */

TEST(uRegionNGSTestHerit, PRINT){
    uRegion uTest("chr1", 100, 119);
   EXPECT_NO_THROW( uTest.print(cout));
}
