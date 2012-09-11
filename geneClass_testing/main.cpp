#include "uFormats.h"
#include "uTags.h"
#include "uRegion.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include "utility/utility.h"
#include <time.h>
#include "gtest.h"

using namespace std;

    #define SOMENUMBER 102343

using namespace testing;
class TestEventListenerProxy : public TestEventListener
{
public:
explicit TestEventListenerProxy(TestEventListener* event_listener);
virtual ~TestEventListenerProxy();

virtual void OnTestProgramStart(const UnitTest& unit_test);
virtual void OnTestIterationStart(const UnitTest& unit_test, int iteration);
virtual void OnEnvironmentsSetUpStart(const UnitTest& unit_test);
virtual void OnEnvironmentsSetUpEnd(const UnitTest& unit_test);
virtual void OnTestCaseStart(const TestCase& test_case);
virtual void OnTestStart(const TestInfo& test_info);
virtual void OnTestPartResult(const TestPartResult& result);
virtual void OnTestEnd(const TestInfo& test_info);
virtual void OnTestCaseEnd(const TestCase& test_case);
virtual void OnEnvironmentsTearDownStart(const UnitTest& unit_test);
virtual void OnEnvironmentsTearDownEnd(const UnitTest& unit_test);
virtual void OnTestIterationEnd(const UnitTest& unit_test, int iteration);
virtual void OnTestProgramEnd(const UnitTest& unit_test);

protected:
TestEventListener* listener;
};


class CaseSummaryAndFailurePrinter : public TestEventListenerProxy
{
public:
explicit CaseSummaryAndFailurePrinter(TestEventListener* default_printer)
    : TestEventListenerProxy(default_printer)
{
}

virtual void OnEnvironmentsTearDownStart(const UnitTest& /*unit_test*/) { }
virtual void OnEnvironmentsSetUpStart(const UnitTest& /*unit_test*/) { }
virtual void OnTestStart(const TestInfo& /*test_info*/) { }

virtual void OnTestEnd(const TestInfo& test_info) {
    if (test_info.result()->Failed())
        listener->OnTestEnd(test_info);
    }
};



class SummaryAndFailurePrinter : public CaseSummaryAndFailurePrinter
{
public:
explicit SummaryAndFailurePrinter(TestEventListener* default_printer)
    : CaseSummaryAndFailurePrinter(default_printer)
{
}

virtual void OnTestCaseStart(const TestCase& /*test_case*/) { }
virtual void OnTestCaseEnd(const TestCase& /*test_case*/) { }
};


class uItem : public uGenericNGS{

    int chr;

};

class ourDerivedClass : public uGenericNGSChrom<uGenericNGS> {

public:
  ourDerivedClass():uGenericNGSChrom(){};
  ourDerivedClass(std::string chrom):uGenericNGSChrom(chrom){};
  ourDerivedClass(std::string chrom, long int size):uGenericNGSChrom(chrom,size){};
};

class ourDerivedExperiment : public uGenericNGSExperiment<ourDerivedClass, uGenericNGS> {

};



TEST(ExperimentDivide, DIVIDECHROMINTONBIN){

    ourDerivedExperiment ourExp;
    ourExp.addSite(uGenericNGS("chr21", 100, 199));
    ourExp.addSite(uGenericNGS("chr22", 100, 199));
    ourExp.addSite(uGenericNGS("chr23", 100, 199));
    ourExp.addSite(uGenericNGS("chr24", 100, 199));
    ourExp.divideItemsIntoNBins(4);

    EXPECT_EQ(ourExp.count(),16);
}
TEST(ExperimentDivide, DIVIDECHROMINTOBINOFSIZE){


    ourDerivedExperiment ourExp;

    ourExp.addSite(uGenericNGS("chr21", 100, 199));
    ourExp.addSite(uGenericNGS("chr22", 100, 199));
    ourExp.addSite(uGenericNGS("chr23", 100, 199));
    ourExp.addSite(uGenericNGS("chr24", 100, 199));
    ourExp.divideItemsIntoBinofSize(50);
    EXPECT_EQ(ourExp.count(),8);

}



/**< Testing our base container FUNCTIONS */
/**< Not that we should never implement this without a derived child, so as such the construction is not accesible */
TEST(uGenericNGSChromTest, CTR){

    EXPECT_NO_THROW(ourDerivedClass uChromTest);
    EXPECT_NO_THROW(ourDerivedClass uChromTest2("chr1"));
    EXPECT_NO_THROW(ourDerivedClass uChromTest3(""));
}
TEST(uGenericNGSChromTest, FUNCTIONS){


    ourDerivedClass uChromEmpty;
    ourDerivedClass uChromChr1("chr1");
    ourDerivedClass uChromNull("");
    ourDerivedClass uChromBig("chrBig", SOMENUMBER);

    EXPECT_EQ(uChromEmpty.getChr() , uChromNull.getChr());

    EXPECT_EQ("chr1", uChromChr1.getChr());
    EXPECT_EQ(0, uChromChr1.count());
    EXPECT_EQ(0, uChromChr1.getChromSize());
    EXPECT_ANY_THROW(uChromChr1.setChromSize(-1));
    EXPECT_NO_THROW(uChromChr1.setChromSize(SOMENUMBER));
    EXPECT_EQ(uChromBig.getChromSize(),uChromChr1.getChromSize());

    uChromBig.setChr("chr3");
    EXPECT_EQ("chr3", uChromBig.getChr());

    EXPECT_ANY_THROW(uChromBig.getSite(0));
}


TEST(factoryTest, uGenericTest){

    EXPECT_ANY_THROW(factory::makeNGSfromTabString("chr2 400 300"));

    uGenericNGS ourTest=factory::makeNGSfromTabString(("chr2 200 300"));

    EXPECT_EQ(ourTest.getChr(), "chr2");
    EXPECT_EQ(ourTest.getStart(), 200);
    EXPECT_EQ(ourTest.getEnd(), 300);
}

TEST(writeTest, uGenecExpTest){

   // ourDerivedExperiment testExp;
   //uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>,uGenericNGS> testExp;
    uTagsExperiment textExp;
    //testExp.addSite(uGenericNGS("chr21", 300, 800));
    textExp.writeAsBedFile(cout);


}

using namespace testing;
int main(int argc, char **argv) {



  ::testing::InitGoogleTest(&argc, argv);
 return  RUN_ALL_TESTS();;
}
