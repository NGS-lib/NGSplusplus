
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


class StandardRegionChr
{
public:
	StandardRegionChr()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
       severalItems.setChr("chr2");
       severalItems.addData(uRegion("chr2",100,200));
       uRegion withSignal("chr2",100,104,StrandDir::REVERSE,0.2f);
       withSignal.setSignal({1,3,6,2,3});
       severalItems.addData(withSignal);
       withSignal.setSignal({0,0,0,0,0});
       severalItems.addData(withSignal);
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
	uRegionChrom severalItems;
	uRegionChrom threeItems;
    uTagsChrom densitySample;
};


TEST(uRegionCHR_AssigmentOperator, NORMAL){
    StandardRegionChr fixedChrom;
	uRegionChrom newChrom=fixedChrom.severalItems;
    auto curItr= newChrom.begin();
    auto Standarditr= fixedChrom.severalItems.begin();

    for(int i=0; i<newChrom.count(); i++){
        EXPECT_TRUE(curItr->isEqual(*Standarditr));
        curItr++;
        Standarditr++;
     }
}
TEST(uRegionCHR_MeasureDensityOverlap, VALID) {
    StandardRegionChr fixedChrom;
  //  fixedChrom.threeItems.sortSites();
  //  fixedChrom.densitySample.sortSites();
    fixedChrom.threeItems.setChromSize(200000);
    fixedChrom.threeItems.measureDensityOverlap(fixedChrom.densitySample);
    auto itr =fixedChrom.threeItems.begin();
  //  std::cout << fixedChrom.densitySample.count() <<std::endl;
    EXPECT_EQ(7, itr->getCount());
    itr++;
    EXPECT_EQ(3, itr->getCount());
    itr++;
    EXPECT_EQ(0, itr->getCount());

   }
TEST(uRegionCHR_GenerateSignal, VALID) {
    StandardRegionChr fixedChrom;
   fixedChrom.threeItems.sortSites();
   fixedChrom.densitySample.sortSites();
   fixedChrom.threeItems.setChromSize(200000);
   fixedChrom.threeItems.generateSignal(fixedChrom.densitySample);
   ASSERT_FALSE(true);
   //fixedChrom.threeItems.writeSignal(std::cout,'\0');
}


TEST(uRegionCHR_WriteSignal, NOCRASH) {
	 StandardRegionChr fixedChrom;
	 EXPECT_NO_THROW(fixedChrom.severalItems.writeSignal(std::cout,'\0'));
}

TEST(uRegionCHRTest_copyCtr, NORMAL){

    StandardRegionChr fixedChrom;
	uRegionChrom newChrom(fixedChrom.severalItems);
    auto curItr= newChrom.begin();
    auto Standarditr= fixedChrom.severalItems.begin();

    for(int i=0; i<newChrom.count(); i++){
        EXPECT_TRUE(curItr->isEqual(*Standarditr));
        curItr++;
        Standarditr++;
     }

}
