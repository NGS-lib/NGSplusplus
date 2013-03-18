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
using namespace NGS;
const string HG18UCSC="../data/GENEPRED/hg18GenePred.genePred";

TEST(uGeneExpTest_addData, ParseNormal)
{
    uGeneExperiment testExp;
    uParser genePredParser(HG18UCSC,"GENEPRED");
    EXPECT_NO_THROW(testExp.loadWithParser(genePredParser,1));
    EXPECT_EQ(1, testExp.count());

    auto chr1P= testExp.getpChrom("chr1");
    EXPECT_EQ(4,chr1P->begin()->featureCount());

}

TEST(uGeneExpTest_addData, PARSEIF)
{
    uGeneExperiment testExp;
    uParser genePredParser(HG18UCSC,"GENEPRED");
    int yaya=0;
    /**< Skip first */
    EXPECT_NO_THROW(testExp.loadWithParser_if(genePredParser,[&](const uGene item){
                                              if (yaya>=1)
                                                return true;
                                              yaya++;
                                            return false;
                                               },  1));
    EXPECT_EQ(1, testExp.count());
    auto chr1P= testExp.getpChrom("chr1");
    EXPECT_EQ(3,chr1P->begin()->featureCount());

}

TEST(uGeneExpTest_addData, VIRTUALTEST)
{
    uParser genePred(HG18UCSC,"GENEPRED");
    uGeneExperiment curexp;

    EXPECT_NO_THROW(curexp.addData(genePred.getNextEntry()));
    EXPECT_EQ(1, curexp.count());
    auto chr1P= curexp.getpChrom("chr1");
    EXPECT_EQ(4,chr1P->begin()->featureCount());

}

