
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



