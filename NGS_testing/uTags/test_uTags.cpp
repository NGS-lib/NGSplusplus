
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

class StandardTags
{
public:
	StandardTags()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
       simpleTag.setChr("chr1");
       simpleTag.setStartEnd(500,800);

       PETag.setChr("chr1");
       PETag.setStartEnd(500,600);
       PETag.setPELength(200);

       CompletTag.setChr("chr1");
       CompletTag.setStartEnd(400,405);
       CompletTag.setSequence("GAGACG");
       CompletTag.setCigar("6M");
	}

    uTags simpleTag;
    uTags PETag;
    uTags CompletTag;
};
TEST(uTagsTest_getCopy, VALID){

   StandardTags myTags;
   auto newTags=myTags.CompletTag.getCopy();
   EXPECT_TRUE( myTags.CompletTag.isEqual(newTags));
}
TEST(uTagsTest_getCompletedCopy, VALID) {

	StandardTags myTags;
	auto completedPE= myTags.PETag.getCompletedCopy();
	EXPECT_EQ( (301),completedPE.getLength());
	EXPECT_EQ( (600+200),completedPE.getEnd());
    EXPECT_EQ(500,completedPE.getStart());
}
TEST(uTagsTest_isEqual, VALIDCOPY) {

    StandardTags myTags;
	EXPECT_TRUE(myTags.CompletTag.isEqual(myTags.CompletTag));

}
TEST(uTagsTest_createToken, VALID) {
    StandardTags myTags;
	ASSERT_NO_THROW(myTags.CompletTag.createToken());
    uToken newToken=myTags.CompletTag.createToken();
	ASSERT_NO_THROW(uTags fromToken(newToken));
	uTags fromToken(newToken);
	EXPECT_TRUE(fromToken.isEqual(myTags.CompletTag));
}

TEST(uTagsTest_copyCtr, NORMAL){

     StandardTags myTags;
    uTags newTag(myTags.CompletTag);
    EXPECT_TRUE(newTag.isEqual(myTags.CompletTag));
}


TEST(uTagsTest_getCompletedTag, FORWARD){
    uTags  PETag("chr1",100,200);
    PETag.setPELength(50);
    EXPECT_EQ(PETag.getCompletedCopy().getEnd(),250);

}
TEST(uTagsTest_getCompletedTag, BACKWARD){
    uTags  PETag("chr1",100,200,StrandDir::REVERSE);
    PETag.setPELength(50);
    EXPECT_EQ(PETag.getCompletedCopy().getStart(),50);

}

TEST(uTagsTest_ctr, POLYREG){
   uRegion simpleBasic("chr1",100,200);
   ASSERT_NO_THROW(uTags fromBasic(simpleBasic) );
   uRegion regFrom("chr1", 100, 200, 0.4f);
   regFrom.setScore(0.3f,2);
   regFrom.setSignal(4, 0.2f);

   uTags fromBasic(regFrom);

   EXPECT_EQ( fromBasic.getChr(),regFrom.getChr() ) ;
   EXPECT_EQ(fromBasic.getStart(),regFrom.getStart());
   EXPECT_EQ(fromBasic.getEnd(),regFrom.getEnd());
   EXPECT_EQ(fromBasic.getScoreVector(),regFrom.getScoreVector());
}

TEST(uTagsTest_ctr, POLYBASIC){

   ASSERT_NO_THROW(uTags fromBasic(uBasicNGS("chr1", 100, 200)) );
   uBasicNGS basicFrom("chr1", 100, 200, 0.4f);
   basicFrom.setScore(0.3f,2);
   uTags fromBasic(basicFrom);

   EXPECT_EQ(fromBasic.getChr(),basicFrom.getChr());
   EXPECT_EQ(fromBasic.getStart(),basicFrom.getStart());
   EXPECT_EQ(fromBasic.getEnd(),basicFrom.getEnd());
   EXPECT_EQ(fromBasic.getScoreVector(),basicFrom.getScoreVector());
}


TEST(uTagsTest_ctr, POLYGENE){

   ASSERT_NO_THROW(uTags fromBasic(uGene("chr1", 100, 200)) );
   uGene basicFrom("chr1", 100, 200, 0.4f);
   basicFrom.setScore(0.3f,2);
   uTags fromBasic(basicFrom);

   EXPECT_EQ(fromBasic.getChr(),basicFrom.getChr());
   EXPECT_EQ(fromBasic.getStart(),basicFrom.getStart());
   EXPECT_EQ(fromBasic.getEnd(),basicFrom.getEnd());
   EXPECT_EQ(fromBasic.getScoreVector(),basicFrom.getScoreVector());
}

TEST(uTagsTest_ctr, DEFAULT){
    uTags uTest;
    EXPECT_EQ(StrandDir::FORWARD,uTest.getStrand());
    EXPECT_FALSE(uTest.isPE());
    EXPECT_EQ(0,uTest.getStart());
    EXPECT_EQ(0,uTest.getEnd());
}

TEST(uTagsTest_ctr, 3OAR){
    uTags Utest("chr1", 100, 200);
    EXPECT_EQ("chr1",Utest.getChr());
    EXPECT_EQ(StrandDir::FORWARD,Utest.getStrand());
    EXPECT_FALSE(Utest.isPE());
    EXPECT_EQ(100,Utest.getStart());
    EXPECT_EQ(200,Utest.getEnd());
}
// TODO: Copy operator
TEST(uTagsTest_ctr, ASSIGNMENT){
    uTags copyFrom("chr1", 100, 200);
    uTags copyTo;
    uTags copyFilled("chr4", 100000, 2000000);
    copyFrom.setName("My Name is");
    copyFrom.setPhred("My Phred is");
    copyFrom.setCigar("Hohoho");
    copyFrom.setPELength(10);
    copyFrom.setMapQual(33);
    copyFrom.setSequence("LalaBonjouLalaBonjouLalaBonjouLalaBonjouLalaBonjouLalaBonjouLalaBonjouLalaBonjouLalaBonjouLalaBonjoud");
    copyTo.setChr("chr22");

    copyTo.setMapped(true);
    copyTo=copyFrom;
    EXPECT_EQ(copyFrom.getPeLength(), copyTo.getPeLength());
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

TEST(uTagsTest_ctr, 3ParCTRMinus){
   EXPECT_ANY_THROW(uTags Utest("chr1", 200, 100));
}

TEST(uTags_getToken, GETTOKENTEST){
    uTags uTest("chr1", 100, 119);
    uTest.setSequence("nnnnnnnnnnnnnnnnnnnn");
    uTest.setCigar("20M");
    uTest.setFlag(132);
    using namespace utility;
    uToken ourToken= uTest.createToken();
    EXPECT_EQ("chr1", ourToken.getParam(token_param::CHR));
    EXPECT_EQ("100", ourToken.getParam(token_param::START_POS));
    EXPECT_EQ("119", ourToken.getParam(token_param::END_POS));
    EXPECT_EQ("+", ourToken.getParam(token_param::STRAND));
    EXPECT_EQ("20M", ourToken.getParam(token_param::CIGAR));
    EXPECT_EQ(132, utility::stoi(ourToken.getParam(token_param::FLAGS)));

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

TEST(uTagsExpTest, distinctTest)
{
    uTagsExperiment emptyExp;
  //  emptyExp.sortSites();
    auto chrCopy=emptyExp.getDistinct("chr1",10, 20);
    EXPECT_EQ(0,chrCopy.count());

    emptyExp.addData(uTags("chr5",800, 1400));
    emptyExp.addData(uTags("chr2",400, 2000));
    emptyExp.addData(uTags("chr5",400, 2200));
    emptyExp.addData(uTags("chr5",1000, 3000));
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
    EXPECT_EQ(ourTest.getLength(),40);
    EXPECT_EQ(ourTest.getStart(),9719905);
    EXPECT_EQ(ourTest.getEnd(),9719944);
    EXPECT_EQ(ourTest.getStrand(),StrandDir::FORWARD);

    uTags minusStrand;
    minusStrand=factory::makeTagfromSamString("HWI-ST333_0111_FC:6:2202:20769:154221#TAGTCG/3	83	chr21	9719905	15	40M	=	9719985	120	AGCAATTATCTTCACATAAAAACTACACAGAAACTTTCTG	aaacceeegggggiiiiiiiiiihiihihiiiiiiiiiih	X0:i:2	X1:i:57	MD:Z:2G2A27T6	XG:i:0	AM:i:0	NM:i:3	SM:i:0	XM:i:3	XO:i:0	XT:A:R");

    EXPECT_EQ(minusStrand.getStrand(),StrandDir::REVERSE);

}
