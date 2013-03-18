
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
   regExp.measureDensityOverlap(tagExp);
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


