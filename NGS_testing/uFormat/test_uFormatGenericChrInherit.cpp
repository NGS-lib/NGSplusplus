/**< Test Common inherited functions */
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
class ChromDivide : public testing::Test {
 protected:
/**< As always, this is inclusive so 100-199 is of size 100 */
  virtual void SetUp() {

    uChromTestOverlap.setChr("chr1");
    uChromTestOverlap.addData(uBasicNGS ("chr1", 300, CHROMDIVIDESIZE));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 199));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 299));
  }
uBasicNGSChrom uChromTestOverlap;
};


/**< Apply on AllVecData */
TEST(uBasicNGSCHR_applyAndGetVecData, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_applyAndGetVecData, ALLMODIFIED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_applyAndGetVecData, ONEMPTY){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_applyAndGetVecData, THROWUGENE){
       ASSERT_TRUE(false);
 }
/**<  computeOnAllSites*/
TEST(uBasicNGSCHR_computeOnAllSites, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_computeOnAllSites, ALLMODIFIED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_computeOnAllSites, ONEMPTY){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSCHR_computeOnAllSites, ONEMPTY){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSCHR_computeOnAllSites, THROWUGENE){
       ASSERT_TRUE(false);
 }

/**<  getSpecificSites*/
TEST(uBasicNGSCHR_getSpecificSites, NONECOUNTED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_getSpecificSites, SOMECOUNTED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_getSpecificSites, ALLCOUNTED){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSCHR_getSpecificSites, EMPTY){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSCHR_getSpecificSites, THROWUGENE){
       ASSERT_TRUE(false);
 }
/**<  removeSpecificSites
 */
TEST(uBasicNGSCHR_removeSpecificSites, NONEREMOVED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_removeSpecificSites, SOMEREMOVED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_removeSpecificSites, ALLREMOVED){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSCHR_removeSpecificSites, EMPTY){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSCHR_removeSpecificSites, THROWUGENE){
       ASSERT_TRUE(false);
 }

/**< applyOnAllSites */
TEST(uBasicNGSCHR_applyOnAllSites, EMPTY){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_applyOnAllSites, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_applyOnAllSites, EXCEPTION){
       ASSERT_TRUE(false);
 }
/**< applyOnAllSites */
TEST(uBasicNGSCHR_applyOnAllSitesConst, EMPTY){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_applyOnAllSitesConst, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_applyOnAllSitesConst, EXCEPTION){
       ASSERT_TRUE(false);
 }
/**< accumulateSitesInfo */
TEST(uBasicNGSCHR_accumulateSitesInfos, EMPTY){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_accumulateSitesInfo, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHR_accumulateSitesInfo, EXCEPTION){
       ASSERT_TRUE(false);
 }
/**< MinMaxSites */
TEST(uBasicNGSCHR_minAndMaxSites, NORMAL){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSCHR_minAndMaxSites, EMPTY){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSCHR_minAndMaxSites, EXCEPTION){
       ASSERT_TRUE(false);
 }

/**<countSitesWithProperty */
TEST(uBasicNGSCHR_countSitesWithProperty, NORMAL){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSCHR_countSitesWithProperty, EMPTY){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSCHR_countSitesWithProperty, EXCEPTION){
       ASSERT_TRUE(false);
 }
/**< applyOnAllChroms */


 TEST(uBasicNGSCHR_applyOnAllChroms, VALID){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSCHR_applyOnAllChroms, EMPTY){
       ASSERT_TRUE(false);
 }
