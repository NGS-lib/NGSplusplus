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

class StandardGenes
{
public:
	StandardGenes()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
       simpleGene.setChr("chr1");
       simpleGene.setStartEnd(500,800);

       simpleGeneID.setChr("chr1");
       simpleGeneID.setStartEnd(500,800);
       simpleGeneID.setID("uasd.23");
       simpleGeneID.setTranscript("uasd.233sw");

       manyFeaturesGene.setChr("chr1");
       manyFeaturesGene.setStartEnd(400,1000);
       manyFeaturesGene.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"GENES","2");
       manyFeaturesGene.addFeature(75, 150,StrandDir::FORWARD,featureType::EXON,"WEIRD","2");
       manyFeaturesGene.addFeature(3000, 4000,StrandDir::FORWARD,featureType::EXON,"COPY","2");
	}
	uGene manyFeaturesGene;
    uGene simpleGene;
    uGene simpleGeneID;
};


TEST(uGenesTest_ctr, POLYREG){
   uBasicNGS simpleBasic("chr1",100,200);
   ASSERT_NO_THROW(uGene fromBasic(simpleBasic) );
   uRegion regFrom("chr1", 100, 200, 0.4f);
   regFrom.setScore(0.3f,2);
   regFrom.setSignal(4, 0.2f);

   uGene fromBasic(regFrom);

   EXPECT_EQ( fromBasic.getChr(),regFrom.getChr() ) ;
   EXPECT_EQ(fromBasic.getStart(),regFrom.getStart());
   EXPECT_EQ(fromBasic.getEnd(),regFrom.getEnd());
   EXPECT_EQ(fromBasic.getScoreVector(),regFrom.getScoreVector());
}

TEST(uGenesTest_ctr, POLYBASIC){

   ASSERT_NO_THROW(uGene fromBasic(uBasicNGS("chr1", 100, 200)) );
   uBasicNGS basicFrom("chr1", 100, 200, 0.4f);
   basicFrom.setScore(0.3f,2);
   uGene fromBasic(basicFrom);

   EXPECT_EQ(fromBasic.getChr(),basicFrom.getChr());
   EXPECT_EQ(fromBasic.getStart(),basicFrom.getStart());
   EXPECT_EQ(fromBasic.getEnd(),basicFrom.getEnd());
   EXPECT_EQ(fromBasic.getScoreVector(),basicFrom.getScoreVector());
}

TEST(uGenesTest_ctr, POLYTAGSs){

   ASSERT_NO_THROW(uGene fromTag(uTags("chr1", 100, 200)) );
   uTags tagFrom("chr1", 100, 200, 0.4f);
   tagFrom.setScore(0.3f,2);
   uGene fromTag(tagFrom);

   EXPECT_EQ(fromTag.getChr(),tagFrom.getChr());
   EXPECT_EQ(fromTag.getStart(),tagFrom.getStart());
   EXPECT_EQ(fromTag.getEnd(),tagFrom.getEnd());
   EXPECT_EQ(fromTag.getScoreVector(),tagFrom.getScoreVector());
}


TEST(uGeneTest_setgetTranscript, VALID){
    StandardGenes myGenes;
    EXPECT_EQ("",myGenes.simpleGene.getTranscript());
    EXPECT_NO_THROW(myGenes.simpleGene.setTranscript("233.w"));
    EXPECT_EQ("233.w",myGenes.simpleGene.getTranscript());
    EXPECT_NO_THROW(myGenes.simpleGene.setTranscript(""));
    EXPECT_EQ("",myGenes.simpleGene.getTranscript());
}

TEST(uGeneTest_setgetID, VALID){
    StandardGenes myGenes;
    EXPECT_NO_THROW(myGenes.simpleGene.setID("JERA2"));
    EXPECT_EQ("JERA2",myGenes.simpleGene.getID());
    EXPECT_NO_THROW(myGenes.simpleGene.setID(""));
    EXPECT_EQ("",myGenes.simpleGene.getID());
}



TEST(uGeneTest_addFeature, SIMPLEFEATURE){
    StandardGenes myGenes;
   EXPECT_NO_THROW( myGenes.manyFeaturesGene.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"",""));
}

TEST(uGeneTest_addFeature, INVALIDFEATURE){
    StandardGenes myGenes;
    EXPECT_ANY_THROW(myGenes.manyFeaturesGene.addFeature(-200, 300,StrandDir::FORWARD,featureType::EXON,"",""));
    EXPECT_ANY_THROW(myGenes.manyFeaturesGene.addFeature(200, 100,StrandDir::FORWARD,featureType::EXON,"",""));
    EXPECT_ANY_THROW(myGenes.manyFeaturesGene.addFeature(-200, -300,StrandDir::FORWARD,featureType::EXON,"",""));
}

TEST(uGeneTest_addFeature, COMPLETEFEATURES){
   StandardGenes myGenes;
   EXPECT_NO_THROW( myGenes.manyFeaturesGene.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"GENES","243"));
   EXPECT_NO_THROW( myGenes.manyFeaturesGene.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"WEIRD","244"));
   EXPECT_NO_THROW( myGenes.manyFeaturesGene.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"COPY","242"));
}


TEST(uGeneTest_isOverlappingFeature, DOESOVERLAP){
   StandardGenes myGenes;
   EXPECT_TRUE( myGenes.manyFeaturesGene.isOverlappingFeature(25, 50));
   EXPECT_TRUE( myGenes.manyFeaturesGene.isOverlappingFeature(300 , 300));
   EXPECT_TRUE( myGenes.manyFeaturesGene.isOverlappingFeature(2000, 6000));

}

TEST(uGeneTest_isOverlappingFeature, BETWEENFEATURES){
   StandardGenes myGenes;
   EXPECT_FALSE( myGenes.manyFeaturesGene.isOverlappingFeature(2000,2500));
}

TEST(uGeneTest_isOverlappingFeature, NOFEATURES){
   StandardGenes myGenes;
   EXPECT_FALSE( myGenes.simpleGene.isOverlappingFeature(0, 10000000));
}

TEST(uGeneTest_isOverlappingFeature, SPECIFICDOESOVERLAP){
   StandardGenes myGenes;
   EXPECT_TRUE( myGenes.manyFeaturesGene.isOverlappingFeature(25, 50, featureType::EXON));
}

TEST(uGeneTest_isOverlappingFeature, SPECIFICDOESNOTOVERLAP){
   StandardGenes myGenes;
   EXPECT_FALSE( myGenes.manyFeaturesGene.isOverlappingFeature(25, 50, featureType::CODING));
}

TEST(uGeneTest_isOverlappingFeature, NOFEATURESSPECIFIC){
   StandardGenes myGenes;
   EXPECT_FALSE( myGenes.simpleGene.isOverlappingFeature(0, 50000, featureType::EXON));
}

TEST(uGeneTest_isEqual, NOFEATURESSPECIFIC){
   StandardGenes myGenes;
   auto newGene= myGenes.simpleGene;
   newGene.addFeature(-0, 300,StrandDir::FORWARD,featureType::EXON,"","");
   EXPECT_FALSE( myGenes.simpleGene.isEqual(newGene));
}
TEST(uGeneTest_isEqual, SAMEMANY){
   StandardGenes myGenes;
   EXPECT_TRUE( myGenes.manyFeaturesGene.isEqual(myGenes.manyFeaturesGene));

}
TEST(uGeneTest_isEqual, DIFFERENTNOTFEATURES){
   StandardGenes myGenes;
   auto newGene=myGenes.manyFeaturesGene;
    newGene.setID("hahaha");
   EXPECT_FALSE( myGenes.simpleGene.isEqual(newGene));
}

TEST(uGeneTest_getCopy, VALID){
   StandardGenes myGenes;
   auto newGene=myGenes.manyFeaturesGene.getCopy();
   EXPECT_TRUE( myGenes.manyFeaturesGene.isEqual(newGene));
}

TEST(uGeneTest_HasFeatureType, HAS){
   StandardGenes myGenes;
   EXPECT_TRUE( myGenes.manyFeaturesGene.hasFeatureType(featureType::EXON));
}
TEST(uGeneTest_HasFeatureType, EMPTY){
   StandardGenes myGenes;
   EXPECT_FALSE( myGenes.simpleGene.hasFeatureType(featureType::EXON));
}
TEST(uGeneTest_HasFeatureType, DOESNOTHAVE){
   StandardGenes myGenes;
   EXPECT_FALSE( myGenes.manyFeaturesGene.hasFeatureType(featureType::CODING));
}

TEST(uGeneTest_removeFeature, REMOVEONE){
    StandardGenes myGenes;
    auto itr = myGenes.manyFeaturesGene.featureBegin();
    itr ++;
    myGenes.manyFeaturesGene.removeFeature(itr);
    EXPECT_EQ(2, myGenes.manyFeaturesGene.featureCount());
}

TEST(uGeneTest_removeFeature, REMOVEALL){
    StandardGenes myGenes;
    auto itrstart = myGenes.manyFeaturesGene.featureBegin();
    auto itrEnd= myGenes.manyFeaturesGene.featureEnd();
    myGenes.manyFeaturesGene.removeFeature(itrstart,itrEnd);
    EXPECT_EQ(0, myGenes.manyFeaturesGene.featureCount());
}

TEST(uGeneTest_createToken, VALID){
    ASSERT_TRUE(false);
}

TEST(uGeneTest_getFeatureInt, VALID){
    ASSERT_TRUE(false);
}

TEST(uGeneTest_Assignment, VALID){
    ASSERT_TRUE(false);
}
