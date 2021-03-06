
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
/**< Testing our unitary Tag */
TEST(uRegionTest, DefaultCTR){
    uRegion uTest;
    EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
    EXPECT_EQ(0,uTest.getStart());
    EXPECT_EQ(0,uTest.getEnd());
}

TEST(uRegionNGSTestHerit_AssigmentOperator, VALID){

	uRegion firsReg("Chr1",100,103);
	firsReg.setSignal({2,4,2,5});
	firsReg.setIdent("hihi");
	uRegion secondReg=firsReg;
    EXPECT_TRUE(firsReg.isEqual(secondReg));
}


TEST(uRegionTest, UsefulCTR){
    uRegion uTest("chr1",100,200);
    EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
    EXPECT_EQ(100,uTest.getStart());
    EXPECT_EQ(200,uTest.getEnd());
}

TEST(uRegionChrTest, DIVIDEINTONBINFAILCHROM){
    uRegionChrom emptyChrom("chr1");
    emptyChrom.addData(uRegion("chr1",100,200));
    EXPECT_ANY_THROW(emptyChrom.divideItemsIntoNBins(7));

}

TEST(uRegionExpTest, MEASUREDENSITY){
   uRegionExperiment regExp;
   uTagsExperiment tagExp;
   uBasicNGSExperiment basicEXP;
   EXPECT_NO_THROW(regExp.measureDensityOverlap(tagExp));
}


TEST(uRegionTest_ctr, POLYGENE){
   uRegion simpleRegion("chr1",100,200);
   ASSERT_NO_THROW(uGene fromGene(simpleRegion) );
   uGene geneFrom("chr1", 100, 200, 0.4f);
   geneFrom.setScore(0.3f,2);

   uRegion fromGene(geneFrom);

   EXPECT_EQ(fromGene.getChr(),geneFrom.getChr() ) ;
   EXPECT_EQ(fromGene.getStart(),geneFrom.getStart());
   EXPECT_EQ(fromGene.getEnd(),geneFrom.getEnd());
   EXPECT_EQ(fromGene.getScoreVector(),geneFrom.getScoreVector());
}

TEST(uRegionTest_ctr, POLYBASIC){

   ASSERT_NO_THROW(uRegion fromBasic(uBasicNGS("chr1", 100, 200)) );
   uBasicNGS basicFrom("chr1", 100, 200, 0.4f);
   basicFrom.setScore(0.3f,2);
   uRegion fromBasic(basicFrom);

   EXPECT_EQ(fromBasic.getChr(),basicFrom.getChr());
   EXPECT_EQ(fromBasic.getStart(),basicFrom.getStart());
   EXPECT_EQ(fromBasic.getEnd(),basicFrom.getEnd());
   EXPECT_EQ(fromBasic.getScoreVector(),basicFrom.getScoreVector());
}


TEST(uRegionTest_ctr, POLYTAGSs){

   ASSERT_NO_THROW(uRegion fromBasic(uTags("chr1", 100, 200)) );
   uTags tagFrom("chr1", 100, 200, 0.4f);
   tagFrom.setScore(0.3f,2);
   uRegion fromTag(tagFrom);

   EXPECT_EQ(fromTag.getChr(),tagFrom.getChr());
   EXPECT_EQ(fromTag.getStart(),tagFrom.getStart());
   EXPECT_EQ(fromTag.getEnd(),tagFrom.getEnd());
   EXPECT_EQ(fromTag.getScoreVector(),tagFrom.getScoreVector());
}


TEST(uRegionTest_isEQUAL, EQUAL){

	uRegion firsReg("Chr1",100,103);
	firsReg.setSignal({2,4,2,5});
	uRegion otherReg("Chr1",100,103);
	otherReg.setSignal({2,4,2,5});

	EXPECT_TRUE(firsReg.isEqual(otherReg));
}


TEST(uRegionTest_isEQUAL, SIGNALDIFF){
	uRegion firsReg("Chr1",100,103);
	firsReg.setSignal({2,4,2,5});
	uRegion otherReg("Chr1",100,103);
	otherReg.setSignal({2,2,2,5});

	EXPECT_FALSE(firsReg.isEqual(otherReg));

}

TEST(uRegionTest_isEQUAL, OTHERDIFF){

	uRegion firsReg("Chr1",100,103);
	firsReg.setSignal({2,4,2,5});
	uRegion otherReg("Chr1",100,103);
	otherReg.setSignal({2,4,2,5});
    otherReg.setDensity(4.2f);

	EXPECT_FALSE(firsReg.isEqual(otherReg));
	uRegion thirdReg=firsReg;
	firsReg.setIdent("yor");
    EXPECT_FALSE(firsReg.isEqual(thirdReg));
}

TEST(uRegionTest_copyCtr, NORMAL){

    uRegion basicRegWithSignal("chr1", 100, 103);
    basicRegWithSignal.setSignal({2,3,6,3});
    uRegion newReg(basicRegWithSignal);

    EXPECT_EQ(newReg.getSignal(), basicRegWithSignal.getSignal());
    ASSERT_TRUE(newReg.isEqual(basicRegWithSignal));

}

TEST(uRegionTest_SetGetSignal, INVALID){

	uRegion emptyReg;
	EXPECT_ANY_THROW(emptyReg.setSignal({2,4,3,5,2}));

}
TEST(uRegionTest_SetGetSignal, VALID){

	uRegion emptyReg("chr1", 100, 104);
	emptyReg.setSignal({2,4,3,5,2});
    EXPECT_EQ(std::vector<float>({2,4,3,5,2}), emptyReg.getSignal());
}

TEST(uRegionTest_WriteSignal, VALID){
	uRegion emptyReg("chr1", 100, 103);
	emptyReg.setSignal({2,4,3,5});
	EXPECT_NO_THROW(emptyReg.writeSignal(std::cerr,'\0'));
}
TEST(uRegionTest_AssigmentOperator, VALID){

    uRegion basicRegWithSignal("chr1", 100, 103);
    basicRegWithSignal.setSignal({2,3,6,3});
    uRegion newReg=basicRegWithSignal;

    EXPECT_EQ(newReg.getSignal(), basicRegWithSignal.getSignal());
    ASSERT_TRUE(newReg.isEqual(basicRegWithSignal));

}

