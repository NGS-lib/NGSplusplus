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



class StandardGeneChr
{
public:
	StandardGeneChr()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
       severalItems.setChr("chr2");
       severalItems.addData(uGene("chr2",100,200));
       uGene oneData("chr2",100,104,StrandDir::REVERSE,0.2f);
       severalItems.addData(oneData);
       severalItems.addData(oneData);
        threeItems.setChr("chr1");
        threeItems.addData(uRegion("chr1", 100, 200));
        threeItems.addData(uRegion("chr1", 150, 250));
        threeItems.addData(uRegion("chr1", 1000, 2500));

        densitySample.setChr("chr1");
        densitySample.setChromSize(2000000);
        densitySample.addData(uTags("chr1", 75, 100));
        densitySample.addData(uTags("chr1", 76, 105));
        densitySample.addData(uTags("chr1", 75, 100));
        densitySample.addData(uTags("chr1", 75, 101));
        densitySample.addData(uTags("chr1", 95, 100));
        densitySample.addData(uTags("chr1", 95, 170));
        densitySample.addData(uTags("chr1", 160, 210));
        densitySample.addData(uTags("chr1", 205, 210));


	}
	uGeneChrom severalItems;
	uGeneChrom threeItems;
    uTagsChrom densitySample;
};



TEST(uGeneChromTest_vectorCtr, VALID) {
    EXPECT_NO_THROW(uGeneChrom newChrom(vector<uGene>({uGene("chr1",100,200,"haha"),uGene("chr1",100,200,"hihi") }) ));
}

TEST(uGeneChromTest_vectorCtr, INVALID) {
    EXPECT_ANY_THROW(uGeneChrom newChrom( vector<uGene>({uGene("chr1",100,200),uGene("chr1",100,200),uGene("chr1",100,200),uGene("chr4",100,200) }) ));
    EXPECT_ANY_THROW(uGeneChrom newChrom2( vector<uGene>({uGene("chr1",100,200,"haha"),uGene("chr1",100,200,"haha") }) ));
}

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

TEST(uGeneChromTest_addData, VALID)
{
    uGeneChrom newChrom("chr1");
    uParser Parser("../data/GENEPRED/hg18GenePred.genePred", "GENEPRED");
    EXPECT_NO_THROW(newChrom.addData( Parser.getNextEntry()));

}
TEST(uGeneChromTest_addData, DUPLICATE)
{
    uGeneChrom newChrom("chr1");
    uGene newGene("chr1", 100, 200,"MyGene","YAY");
    newChrom.addData(newGene);
    EXPECT_ANY_THROW(newChrom.addData(newGene));
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
const string GTFUCSC="../data/GTF/knownGene.gtf";
const string GTFsmall="../data/GTF/endl.gtf";

TEST(uGeneChromTest_addDataGTF, ParseOne)
{
    uGeneChrom newChrom("chr1");
    uParser Parser(GTFsmall, "GTF");

    ASSERT_NO_THROW(newChrom.addData( Parser.getNextEntry()));


    std::vector<uGene>::const_iterator itr= newChrom.findGene("uc001aaa.3");
    EXPECT_EQ(0,itr->featureCount());
    EXPECT_EQ(11874,itr->getStart());

}
TEST(uGeneChromTest_addDataGTF, PARSE_SEVERAL)
{

    uGeneChrom newChrom("chr1");
    uParser Parser(GTFUCSC, "GTF");
    int count=0;
    while (count!=50){
        ASSERT_NO_THROW(newChrom.addData(Parser.getNextEntry()));
        count++;
    }

    std::vector<uGene>::const_iterator itr= newChrom.findGene("uc001aaa.3");
    EXPECT_EQ(2,itr->featureCount());
    EXPECT_EQ(11874,itr->getStart());

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


    StandardGeneChr fixedChrom;
	uGeneChrom newChrom(fixedChrom.severalItems);
    auto curItr= newChrom.begin();
    auto Standarditr= fixedChrom.severalItems.begin();

    for(int i=0; i<newChrom.count(); i++){
        EXPECT_TRUE(curItr->isEqual(*Standarditr));
        curItr++;
        Standarditr++;
     }


    }

TEST(uGeneChromTest_Assignment, VALID){
    StandardGeneChr fixedChrom;
	uGeneChrom newChrom=fixedChrom.severalItems;
    auto curItr= newChrom.begin();
    auto Standarditr= fixedChrom.severalItems.begin();

    for(int i=0; i<newChrom.count(); i++){
        EXPECT_TRUE(curItr->isEqual(*Standarditr));
        curItr++;
        Standarditr++;
     }

}

TEST(uGeneChromTes_FindGeneTwoString, VALID){
     uGeneChrom newChrom("chr1");
    uGene newGene("chr1", 100, 200,"MyGene","YAY");
    newChrom.addData(newGene);
    newChrom.addData(uGene("chr1", 200, 200,"MyGene3","yoyo"));
    EXPECT_NO_THROW(newChrom.addData(uGene(uGene("chr1", 200, 200,"MyGene3","yarggle"))));

   auto itr= newChrom.findGene("MyGene3","yoyo");
   EXPECT_TRUE(itr->isEqual(uGene("chr1", 200, 200,"MyGene3","yoyo")));

}
/**< Validity tested via assignement operator */
TEST(uGeneChromTes_GetCopy, VALID){
    uGeneChrom newChrom("chr1");
    uGene newGene("chr1", 100, 200,"MyGene","YAY");
    newChrom.addData(newGene);
    EXPECT_NO_THROW(newChrom.getCopy());
}

TEST(uGeneChromTes_FindNext, VALID){
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
       uGeneChrom newChrom("chr1");
       newChrom.addData(manyFeatures);
       newChrom.addData(secondFeatures);
       newChrom.addData(thirdFeatures);
       newChrom.sortSites();
       EXPECT_EQ(newChrom.end(), newChrom.findNextWithFeature(100,featureType::LOOP_END));

      try {

       EXPECT_EQ(400, newChrom.findNextWithFeature(100,featureType::INTRON)->getStart()) ;
       EXPECT_EQ(450, newChrom.findNextWithFeature(410,featureType::EXON)->getStart()) ;
       EXPECT_EQ(3500, newChrom.findNextWithFeature(410,featureType::ENHANCER)->getStart()) ;

    }
    catch(...)
    {
        std::cerr<<"We threw";
    }
    }
//
//TEST(uGeneChromTes_FindPreceding, VALID){
//       /**< Move to class Set-up */
//       uGene manyFeatures;
//       manyFeatures.setChr("chr1");
//       manyFeatures.setStartEnd(400,600);
//       manyFeatures.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"2");
//       manyFeatures.addFeature(75, 150,StrandDir::FORWARD,featureType::EXON,"2");
//       manyFeatures.addFeature(3000, 4000,StrandDir::FORWARD,featureType::PROMOTER,"2");
//
//       uGene secondFeatures;
//       secondFeatures.setChr("chr1");
//       secondFeatures.setStartEnd(450,1000);
//       secondFeatures.addFeature(50, 300,StrandDir::FORWARD,featureType::EXON,"2");
//       secondFeatures.addFeature(75, 150,StrandDir::FORWARD,featureType::EXON,"2");
//        uGene thirdFeatures;
//       thirdFeatures.setChr("chr1");
//       thirdFeatures.setStartEnd(3500,4000);
//       thirdFeatures.addFeature(75, 150,StrandDir::FORWARD,featureType::ENHANCER,"2");
//       uGeneChrom newChrom("chr1");
//       newChrom.addData(manyFeatures);
//       newChrom.addData(secondFeatures);
//       newChrom.addData(thirdFeatures);
//       newChrom.sortSites();
//       EXPECT_EQ(newChrom.end(), newChrom.findPrecedingWithFeature(500000,featureType::LOOP_END));
//
//       EXPECT_EQ(400, newChrom.findPrecedingWithFeature(480,featureType::PROMOTER)->getStart()) ;
//       EXPECT_EQ(450, newChrom.findPrecedingWithFeature(480,featureType::EXON)->getStart()) ;
//
//
//       EXPECT_EQ(3500, newChrom.findPrecedingWithFeature(5000,featureType::ENHANCER)->getStart()) ;
//       EXPECT_EQ(450, newChrom.findPrecedingWithFeature(5000,featureType::EXON)->getStart()) ;
//
//    }


TEST(uGeneChromTes_GetIdCountWithOneOrTwoString, VALID){

    uGeneChrom  newChrom("chr1");
    newChrom.addData(uGene("chr1",100,200,"hihi","hoho"));
    newChrom.addData(uGene("chr1",100,200,"hihi","harhar"));
    newChrom.addData(uGene("chr1",100,200,"hihi","lolo"));
    newChrom.addData(uGene("chr1",100,200,"yarg","lolo"));


    EXPECT_EQ(3,newChrom.getIDCount("hihi"));
    EXPECT_EQ(1,newChrom.getIDCount("hihi","hoho"));
    EXPECT_EQ(1,newChrom.getIDCount("yarg","lolo"));
    EXPECT_EQ(1,newChrom.getIDCount("yarg"));
}


