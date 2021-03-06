
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
using namespace BamTools;
const string BAMREAL="../data/BAM/H2AZ.bam";
vector<int> startPosBam2={56262,61780,64367,67342,71763};

TEST(uTagsExpTest_copyCtr, NORMAL){
    ASSERT_TRUE(false);
}


TEST(uTagsExperiment_loadWithBamReader, COMPAREBAMLOADINGTOPARSER)
{
    uTagsExperiment bamLoad, parserLoad;
    BamReader curReader;
    uParser curParser(BAMREAL,"BAM");
    curReader.Open(BAMREAL);
    EXPECT_NO_THROW(bamLoad.loadWithBamTools(curReader, 10));
    EXPECT_NO_THROW(parserLoad.loadWithParser(curParser,10));

    auto chr10PBam= bamLoad.getpChrom("chr10");
    auto chr10PParser= parserLoad.getpChrom("chr10");
    auto itrBam= chr10PBam->begin();
    auto itrParser= chr10PParser->begin();
    for(size_t i=0;i<startPosBam2.size(); i++)
    {
        EXPECT_EQ(itrBam->getStart(),itrParser->getStart());
        itrBam++;
        itrParser++;
    }

}
TEST(uTagsExperiment_loadWithBamReader, VALIDATESTARTLENGHT)
{
    uTagsExperiment bamLoad, parserLoad;
    BamReader curReader;
    curReader.Open(BAMREAL);
    EXPECT_NO_THROW(bamLoad.loadWithBamTools(curReader, 10));
    EXPECT_EQ(10,bamLoad.count());

    auto chr10P= bamLoad.getpChrom("chr10");
    auto itr= chr10P->begin();
    for(size_t i=0;i<startPosBam2.size(); i++)
    {
        EXPECT_EQ(startPosBam2.at(i),itr->getStart());
        EXPECT_EQ(40,itr->getLength());
        itr++;
    }

}

TEST(uTagsExperiment_getRegionSignal, ONLYBASIC) {
     uTagsExperiment oneExp;
     oneExp.addData(uTags("chr1",1010,1015));
     oneExp.addData(uTags("chr1",102,105));
     EXPECT_EQ(std::vector<float>({0,0,1,1,1,1}), oneExp.getRegionSignal("chr1",100,105,true));
}
