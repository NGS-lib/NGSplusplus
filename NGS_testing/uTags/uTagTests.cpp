
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
TEST(uTagsTest, DefaultCTR){
    uTags uTest;
    EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
    EXPECT_FALSE(uTest.isPE());
    EXPECT_EQ(0,uTest.getStart());
    EXPECT_EQ(0,uTest.getEnd());
}

TEST(uTagsTest, 3ParCTR){
    uTags Utest("chr1", 100, 200);
    EXPECT_EQ("chr1",Utest.getChr());
    EXPECT_EQ(StrandDir::FORWARD,Utest.getStrand());
    EXPECT_FALSE(Utest.isPE());
    EXPECT_EQ(100,Utest.getStart());
    EXPECT_EQ(200,Utest.getEnd());
}

TEST(uTagsTest, 3ParCTRMinus){
   EXPECT_ANY_THROW(uTags Utest("chr1", 200, 100));
}

TEST(uTagsTest, SetGet){
    uTags Utest("chr1", 100, 200);
    uTags copyTag;

    Utest.setName("My Name is");
    EXPECT_EQ("My Name is", Utest.getName());

    EXPECT_EQ("", copyTag.getName());

    Utest.setPhred("My Phred is");
    EXPECT_EQ("My Phred is", Utest.getPhred());

    Utest.setPhred("Now it is yar ar ar");
    EXPECT_EQ("Now it is yar ar ar", Utest.getPhred());

    copyTag=Utest;
    EXPECT_EQ("Now it is yar ar ar", Utest.getPhred());

    Utest.setCigar("Hohoho");
    EXPECT_EQ("Hohoho", Utest.getCigar());
}

TEST(uTagsTest, CopyCTR){
    uTags copyFrom("chr1", 100, 200);
    uTags copyTo;
    uTags copyFilled("chr4", 100000, 2000000);


    copyFrom.setName("My Name is");
    copyFrom.setPhred("My Phred is");
    copyFrom.setCigar("Hohoho");
    copyFrom.setPELenght(10);
    copyFrom.setMapQual(33);
    copyFrom.setSequence("LalaBonjourPatate");
    copyTo.setChr("chr22");

    copyTo.setMapped(true);

    copyTo=copyFrom;
    EXPECT_EQ(copyFrom.getPeLenght(), copyTo.getPeLenght());
    EXPECT_EQ(copyFrom.getSequence(), copyTo.getSequence());
    EXPECT_EQ(copyFrom.getMapQual(), copyTo.getMapQual());
    EXPECT_EQ(copyFrom.getStrand(), copyTo.getStrand());
    EXPECT_EQ(copyFrom.getName(), copyTo.getName());
    EXPECT_EQ(copyFrom.getPhred(), copyTo.getPhred());
    EXPECT_EQ(copyFrom.getSequence(), copyTo.getSequence());
    EXPECT_EQ(copyFrom.getStart(), copyTo.getStart());
    EXPECT_EQ(copyFrom.getChr(), copyTo.getChr());
    EXPECT_EQ(copyFrom.getEnd(), copyTo.getEnd());
    EXPECT_EQ(copyFrom.getCigar(), copyTo.getCigar());
    copyTo.setChr("chr21");
    EXPECT_NE(copyTo.getChr(),copyFrom.getChr());
    copyFilled = copyFrom;
    EXPECT_EQ(copyFrom.getName(), copyFilled.getName());
    EXPECT_EQ(copyFrom.getPhred(), copyFilled.getPhred());
    EXPECT_EQ(copyFrom.getSequence(), copyFilled.getSequence());
    EXPECT_EQ(copyFrom.getStart(), copyFilled.getStart());
    EXPECT_EQ(copyFrom.getChr(), copyFilled.getChr());
    EXPECT_EQ(copyFrom.getEnd(), copyFilled.getEnd());
}



/**< Testing DivideItemsIntoNBin */
TEST(uTagsChrTest, DIVIDEINTONBINFAILCHROM){
    uTagsChrom emptyChrom("chr1");
    emptyChrom.addSite(uTags("chr1",100,200));
    EXPECT_ANY_THROW(emptyChrom.divideItemsIntoNBins(7));
    //EXPECT_ANY_THROW(uChromTestOverlap.divideItemsIntoNBins(32));
}



TEST(uTagsChrTest, distinctChrTest)
{
    uTagsChrom emptyChrom("chr1");
    emptyChrom.sortSites();
    auto chrCopy=emptyChrom.getDistinct(10, 20);
    EXPECT_EQ(chrCopy.count(), 0);

    emptyChrom.addSite(uTags("chr1",100,200));
    emptyChrom.addSite(uTags("chr1",300,600));
    emptyChrom.addSite(uTags("chr1",300,700));
    emptyChrom.sortSites();
    std::cerr << emptyChrom.count() <<std::endl;
    chrCopy=emptyChrom.getDistinct(295, 500);
    EXPECT_EQ(1,chrCopy.count());

    emptyChrom.addSite(uTags("chr1",300,600));
    emptyChrom.addSite(uTags("chr1",300,700));
    emptyChrom.sortSites();
    chrCopy=emptyChrom.getDistinct(0, 1);
    EXPECT_EQ(5,chrCopy.count());


    uTagsChrom newChrom("chr5");
    newChrom.addSite(uTags("chr5",800, 1400));
    newChrom.addSite(uTags("chr5",400, 2000));
    newChrom.addSite(uTags("chr5",400, 2200));
    newChrom.addSite(uTags("chr5",1000, 3000));
    newChrom.sortSites();
    uTagsChrom yar;
    yar=newChrom.getDistinct(300, 900);
     EXPECT_EQ(1,yar.count());

}


TEST(uTagsExpTest, distinctTest)
{
    uTagsExperiment emptyExp;
  //  emptyExp.sortSites();
    auto chrCopy=emptyExp.getDistinct("chr1",10, 20);
    EXPECT_EQ(0,chrCopy.count());

    emptyExp.addSite(uTags("chr5",800, 1400));
    emptyExp.addSite(uTags("chr2",400, 2000));
    emptyExp.addSite(uTags("chr5",400, 2200));
    emptyExp.addSite(uTags("chr5",1000, 3000));
    emptyExp.sortSites();
    chrCopy=emptyExp.getDistinct("chr5",300, 900);
    EXPECT_EQ(2,chrCopy.count());
    uTagsExperiment newExp;
    newExp=(emptyExp.getDistinct("chr10",300, 900));
    EXPECT_EQ(4,newExp.count());
    newExp=(emptyExp.getDistinct("chr2",100, 800));
    EXPECT_EQ(3,newExp.count());

}

TEST(uTagsExpTest, derivedSubsetTestEmpty){
    uTagsExperiment emptyExp;
    uTagsChrom chromToCopyTo;
    chromToCopyTo = emptyExp.getSubset("chr1", 10 , 10000);

    EXPECT_EQ(chromToCopyTo.count(), 0);
}


TEST(factoryTest, uTagsTest){

    EXPECT_NO_THROW(uTags ourTest=factory::makeTagfromSamString("HWI-ST333_0111_FC:6:2202:20769:154221#TAGTCG/3	163	chr21	9719905	15	40M	=	9719985	120	AGCAATTATCTTCACATAAAAACTACACAGAAACTTTCTG	aaacceeegggggiiiiiiiiiihiihihiiiiiiiiiih	X0:i:2	X1:i:57	MD:Z:2G2A27T6	XG:i:0	AM:i:0	NM:i:3	SM:i:0	XM:i:3	XO:i:0	XT:A:R"));

    uTags ourTest=factory::makeTagfromSamString("HWI-ST333_0111_FC:6:2202:20769:154221#TAGTCG/3	163	chr21	9719905	15	40M	=	9719985	120	AGCAATTATCTTCACATAAAAACTACACAGAAACTTTCTG	aaacceeegggggiiiiiiiiiihiihihiiiiiiiiiih	X0:i:2	X1:i:57	MD:Z:2G2A27T6	XG:i:0	AM:i:0	NM:i:3	SM:i:0	XM:i:3	XO:i:0	XT:A:R");
    EXPECT_EQ(ourTest.getChr(),"chr21");
    EXPECT_EQ(ourTest.getLenght(),40);
    EXPECT_EQ(ourTest.getStart(),9719905);
    EXPECT_EQ(ourTest.getEnd(),9719944);
    EXPECT_EQ(ourTest.getStrand(),StrandDir::FORWARD);

    uTags minusStrand;
    minusStrand=factory::makeTagfromSamString("HWI-ST333_0111_FC:6:2202:20769:154221#TAGTCG/3	83	chr21	9719905	15	40M	=	9719985	120	AGCAATTATCTTCACATAAAAACTACACAGAAACTTTCTG	aaacceeegggggiiiiiiiiiihiihihiiiiiiiiiih	X0:i:2	X1:i:57	MD:Z:2G2A27T6	XG:i:0	AM:i:0	NM:i:3	SM:i:0	XM:i:3	XO:i:0	XT:A:R");

    EXPECT_EQ(minusStrand.getStrand(),StrandDir::REVERSE);

}
