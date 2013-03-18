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

//class StandardGeneChroms
//{
//public:
//	StandardGeneChroms()
//	{
//		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
//
//       manyElemChrom.setChr("chr1");
//
//       uGene aGene("chr1");
//       aGene.setStartEnd(400,1000);
//       aGene.addFeature(50, 300,featureType::EXON,"GENES","2");
//       aGene.addFeature(75, 150,featureType::EXON,"WEIRD","2");
//       aGene.addFeature(3000, 4000,featureType::EXON,"COPY","2");
//
//       manyElemChrom.addData(aGene);
//
//	}
//	uGeneChrom manyElemChrom;
// //  uGene simpleGene;
// //   uGene simpleGeneID;
//};


//findNextGeneWithFeature



TEST(uGeneChromTest_addData, ParseNormal)
{
    uGeneChrom newChrom("chr1");
    uParser Parser("../data/BED/test.bed", "BED");
    EXPECT_NO_THROW(newChrom.addData( Parser.getNextEntry()));

}

TEST(uGeneChromTest_addData, ParseOne)
{
    uGeneChrom newChrom("chr1");
    uParser Parser("../data/GENEPRED/hg18GenePred.genePred", "GENEPRED");
    EXPECT_NO_THROW(newChrom.addData( Parser.getNextEntry()));

}

const string UCSCSAMPLE="../data/GFF/UCSCSample.gff";

TEST(uGeneChromTest_addDataGFF, ParseOne)
{
    uGeneChrom newChrom("chr22");
    uParser Parser(UCSCSAMPLE, "UCSCGFF");
    while (Parser.eof()==false){
        EXPECT_NO_THROW(newChrom.addData( Parser.getNextEntry()));
    }

    std::vector<uGene>::const_iterator itr= newChrom.findGene("touch1");
    EXPECT_EQ(1,itr->featureCount());

}
TEST(uGeneChromTest_addDataGFF, ParseALL)
{
    uGeneChrom newChrom("chr22");
    uParser Parser(UCSCSAMPLE, "UCSCGFF");
    while (Parser.eof()==false){
        EXPECT_NO_THROW(newChrom.addData( Parser.getNextEntry()));
    }

    std::vector<uGene>::const_iterator itr= newChrom.findGene("touch1");
    EXPECT_EQ(1,itr->featureCount());

    itr= newChrom.findGene("touch2");
    EXPECT_EQ(0,itr->featureCount());

}


TEST(uGeneChromTest_addDataGTF, ParseOne)
{
     ASSERT_TRUE(false);
    uGeneChrom newChrom("chr22");
    uParser Parser(UCSCSAMPLE, "UCSCGFF");
    while (Parser.eof()==false){
        EXPECT_NO_THROW(newChrom.addData( Parser.getNextEntry()));
    }

    std::vector<uGene>::const_iterator itr= newChrom.findGene("touch1");
    EXPECT_EQ(1,itr->featureCount());

}
TEST(uGeneChromTest_addDataGTF, ParseALL)
{
    ASSERT_TRUE(false);
    uGeneChrom newChrom("chr22");
    uParser Parser(UCSCSAMPLE, "UCSCGFF");
    while (Parser.eof()==false){
        EXPECT_NO_THROW(newChrom.addData( Parser.getNextEntry()));
    }

    std::vector<uGene>::const_iterator itr= newChrom.findGene("touch1");
    EXPECT_EQ(1,itr->featureCount());

    itr= newChrom.findGene("touch2");
    EXPECT_EQ(0,itr->featureCount());

}



TEST(uGeneChromTest_addDataGENEPRED, ParseAndCheck)
{
    uGeneChrom newChrom("chr1");

    vector<int> values ={3,2,4,4,2,3,9,10,10,6,10,7,10,11};
    vector<int> start={1115,1115,2475,3083};
    /**< +1 for coding start */
    auto increase =[](int i){return (i+1);};
    std::transform(std::begin(values), std::end(values),std::begin(values),increase);
    uParser Parser("../data/GENEPRED/hg18GenePred.genePred", "GENEPRED");
    for( int i=0;i<14; i++)
        newChrom.addData( Parser.getNextEntry());
    auto itr= newChrom.begin();
    for (int i=0; i<14;i++){
        EXPECT_EQ(values.at(i),itr->featureCount());
        itr++;
        }
    itr = newChrom.begin();
    int i=0;
    for(auto featureItr=itr->featureBegin(); featureItr!=itr->featureEnd(); featureItr++){
        EXPECT_EQ(start.at(i),featureItr->getStart());
        i++;

    }

}


TEST(uGeneChromTest_copyCtr, NORMAL){
    ASSERT_TRUE(false);
}


