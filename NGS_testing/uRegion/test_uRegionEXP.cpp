

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
class StandardRegionExp
{
public:
	StandardRegionExp()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */

       severalItems.addData(uRegion("chr1",100,200));
       severalItems.addData(uRegion("chr1",230,300));
       severalItems.addData(uRegion("chr1",120,250));

        uRegion withSignal("chr2",100,104,StrandDir::REVERSE,0.2f);
        withSignal.setSignal({1,3,6,2,3});
        uRegion withoutSignal("chr2",200,205,StrandDir::REVERSE,0.5f);
        severalItems.addData(withSignal);
         withSignal.setSignal({3,0,0,0,0});
        severalItems.addData(withoutSignal);
          severalItems.addData(withSignal);
	}
	uRegionExperiment severalItems;

};

TEST(uRegionExpTest_copyCtr, NORMAL){

    StandardRegionExp newGroup;
    uRegionExperiment newExp(newGroup.severalItems);

    EXPECT_EQ(newGroup.severalItems.count(),newExp.count());
    auto curItr=newExp.begin()->second.begin();
    auto otherItr=newGroup.severalItems.begin()->second.begin();
    EXPECT_TRUE(curItr->isEqual(*otherItr));

//    ASSERT_TRUE(false);
}

TEST(uRegionExpTest_WriteSignal, VALID){

     ostringstream* osTest = new ostringstream(ostringstream::out);
	 StandardRegionExp newGroup;
	 EXPECT_NO_THROW(newGroup.severalItems.writeSignal(cout,'\0'));


}
TEST(uRegionExpTest_AssigmentOperator, ASDF){
	ASSERT_TRUE(false);
}
TEST(uRegionExpTest_MeasureDensityOverlap, ASDF) {
	ASSERT_TRUE(false);
}
TEST(uRegionExpTest_GenerateSignal, ASDF) {
	ASSERT_TRUE(false);
}
