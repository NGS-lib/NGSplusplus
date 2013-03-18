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

#define CHROMDIVIDESIZE 500
using namespace std;
using namespace NGS;

class StandardMultipleChroms
{
public:
	StandardMultipleChroms()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
       basicOneChr.setChr("chr1");
       basicOneChr.addData(uBasicNGS("chr1",100,200));
       tagsOneChr.setChr("chr1");
       tagsOneChr.addData(uTags("chr1",100,200));
       regionOneChr.setChr("chr1");
       regionOneChr.addData(uRegion("chr1",100,200));

       basicManyChr.setChr("chr1");
       basicManyChr.addData(uBasicNGS("chr1",100,200));
       basicManyChr.addData(uBasicNGS("chr1",230,300));
       basicManyChr.addData(uBasicNGS("chr1",120,250));

       tagsManyChr.setChr("chr1");
       tagsManyChr.addData(uTags("chr1",100,200));
       tagsManyChr.addData(uTags("chr1",230,300));
       tagsManyChr.addData(uTags("chr1",120,250));

       regionManyChr.setChr("chr1");
       regionManyChr.addData(uRegion("chr1",100,200));
       regionManyChr.addData(uRegion("chr1",230,300));
       regionManyChr.addData(uRegion("chr1",120,250));


       basicEmptyChr.setChr("chr1");
       tagsEmptyChr.setChr("chr1");
       regionEmptyChr.setChr("chr1");
	}
	uBasicNGSChrom basicOneChr;
    uBasicNGSChrom basicManyChr;
    uBasicNGSChrom basicEmptyChr;
    uTagsChrom tagsOneChr;
    uTagsChrom tagsManyChr;
    uTagsChrom tagsEmptyChr;
    uRegionChrom regionOneChr;
    uRegionChrom regionManyChr;
    uRegionChrom regionEmptyChr;
};

/**< Test Constructors */
TEST(uTagsTestChrom_ctr, POLYREG){
       StandardMultipleChroms ourChroms;
       EXPECT_NO_THROW(uTagsChrom(ourChroms.regionManyChr));
       auto polyCreated=uTagsChrom(ourChroms.regionManyChr);

//        auto itr = polyCreated.begin();
       EXPECT_EQ(polyCreated.getSite(0).getStart(),ourChroms.regionManyChr.getSite(0).getStart());
       EXPECT_EQ(polyCreated.getSite(1).getStart(),ourChroms.regionManyChr.getSite(1).getStart());
       EXPECT_EQ(polyCreated.getSite(1).getChr(),ourChroms.regionManyChr.getSite(1).getChr());
       EXPECT_EQ(polyCreated.getSite(2).getEnd(),ourChroms.regionManyChr.getSite(2).getEnd());

}

TEST(uTagsTestChrom_ctr, POLYBASIC){
       StandardMultipleChroms ourChroms;
       EXPECT_NO_THROW(uBasicNGSChrom(ourChroms.basicManyChr));
       auto polyCreated=uBasicNGSChrom(ourChroms.basicManyChr);
       EXPECT_EQ(polyCreated.getSite(0).getStart(),ourChroms.basicManyChr.getSite(0).getStart());
       EXPECT_EQ(polyCreated.getSite(1).getStart(),ourChroms.basicManyChr.getSite(1).getStart());
       EXPECT_EQ(polyCreated.getSite(1).getChr(),ourChroms.basicManyChr.getSite(1).getChr());
       EXPECT_EQ(polyCreated.getSite(2).getEnd(),ourChroms.basicManyChr.getSite(2).getEnd());

}

TEST(uTagsChrTest_copyCtr, NORMAL){
    ASSERT_TRUE(false);
}


/**< Testing DivideItemsIntoNBin */
TEST(uTagsTestChrom_divideIntoNBin, EMPTY){
    uTagsChrom emptyChrom("chr1");
    emptyChrom.addData(uTags("chr1",100,200));
    EXPECT_ANY_THROW(emptyChrom.divideItemsIntoNBins(7));
}

TEST(uTagsTestChrom_divideIntoNBin, DISTINCTTEST)
{
    uTagsChrom emptyChrom("chr1");
    emptyChrom.sortSites();
    auto chrCopy=emptyChrom.getDistinct(10, 20);
    EXPECT_EQ(chrCopy.count(), 0);

    emptyChrom.addData(uTags("chr1",100,200));
    emptyChrom.addData(uTags("chr1",300,600));
    emptyChrom.addData(uTags("chr1",300,700));
    emptyChrom.sortSites();
    std::cerr << emptyChrom.count() <<std::endl;
    chrCopy=emptyChrom.getDistinct(295, 500);
    EXPECT_EQ(1,chrCopy.count());

    emptyChrom.addData(uTags("chr1",300,600));
    emptyChrom.addData(uTags("chr1",300,700));
    emptyChrom.sortSites();
    chrCopy=emptyChrom.getDistinct(0, 1);
    EXPECT_EQ(5,chrCopy.count());

    uTagsChrom newChrom("chr5");
    newChrom.addData(uTags("chr5",800, 1400));
    newChrom.addData(uTags("chr5",400, 2000));
    newChrom.addData(uTags("chr5",400, 2200));
    newChrom.addData(uTags("chr5",1000, 3000));
    newChrom.sortSites();
    uTagsChrom yar;
    yar=newChrom.getDistinct(300, 900);
     EXPECT_EQ(1,yar.count());

}
