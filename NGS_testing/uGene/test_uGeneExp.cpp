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
    EXPECT_EQ(4,(int)chr1P->begin()->featureCount());

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
    EXPECT_EQ(3,(int)chr1P->begin()->featureCount());

}

TEST(uGeneExpTest_addData, VIRTUALTEST)
{
    uParser genePred(HG18UCSC,"GENEPRED");
    uGeneExperiment curexp;

    EXPECT_NO_THROW(curexp.addData(genePred.getNextEntry()));
    EXPECT_EQ(1, curexp.count());
    auto chr1P= curexp.getpChrom("chr1");
    EXPECT_EQ(4,(int)chr1P->begin()->featureCount());

}

TEST(uGeneExpTest_copyCtr, NORMAL){
    ASSERT_TRUE(false);
}
TEST(uGeneEXPTest_Assignment, VALID){
    ASSERT_TRUE(false);
}

TEST(uGeneEXPTest_GetCopy, VALID){
    ASSERT_TRUE(false);
}

TEST(uGeneEXPTest_FindNextWithFeature, VALID){

        uGene manyFeatures;
       manyFeatures.setChr("chr1");
       manyFeatures.setStartEnd(400,600);
       manyFeatures.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"2");
       manyFeatures.addFeature(75, 150,StrandDir::FORWARD,featureType::EXON,"2");
       manyFeatures.addFeature(3000, 4000,StrandDir::FORWARD,featureType::INTRON,"2");

       uGene secondFeatures;
       secondFeatures.setChr("chr1");
       secondFeatures.setStartEnd(450,1000);
       secondFeatures.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"2");
       secondFeatures.addFeature(75, 150,StrandDir::FORWARD,featureType::INTRON,"2");
        uGene thirdFeatures;
       thirdFeatures.setChr("chr1");
       thirdFeatures.setStartEnd(3500,10000);
       thirdFeatures.addFeature(75, 150,StrandDir::FORWARD,featureType::ENHANCER,"2");
       uGeneExperiment newExp;
       newExp.addData(manyFeatures);
       newExp.addData(secondFeatures);
       newExp.addData(thirdFeatures);
       newExp.sortSites();
       auto chrItr= newExp.getpChrom("chr1");
       EXPECT_EQ(chrItr->end(), newExp.findNextWithFeature("chr1",100,featureType::LOOP_END));
       EXPECT_EQ(400, newExp.findNextWithFeature("chr1",100,featureType::INTRON)->getStart()) ;
       EXPECT_EQ(450, newExp.findNextWithFeature("chr1",410,featureType::EXON)->getStart()) ;
       EXPECT_EQ(3500, newExp.findNextWithFeature("chr1",410,featureType::ENHANCER)->getStart()) ;


}
TEST(uGeneEXPTest_FindNextWithFeature, INVALID){
        uGeneExperiment someEmptyExp;
        EXPECT_ANY_THROW(someEmptyExp.findNextWithFeature("chr4", 5, featureType::ENHANCER));
}


