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
    uChromTestOverlap.addData(uBasicNGS ("chr1", 300, 500));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 199));
    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 299));
  }
uBasicNGSChrom uChromTestOverlap;
};


/**< Apply on AllVecData */
TEST(uBasicNGSGENCHR_applyAndGetVecData, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_applyAndGetVecData, ALLMODIFIED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_applyAndGetVecData, ONEMPTY){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_applyAndGetVecData, THROWUGENE){
       ASSERT_TRUE(false);
 }
/**<  computeOnAllSites*/
TEST(uBasicNGSGENCHR_computeOnAllSites, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_computeOnAllSites, ALLMODIFIED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_computeOnAllSites, ONEMPTY){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSGENCHR_computeOnAllSites, THROWUGENE){
       ASSERT_TRUE(false);
 }

/**<  getSpecificSites*/
TEST(uBasicNGSGENCHR_getSpecificSites, NONECOUNTED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_getSpecificSites, SOMECOUNTED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_getSpecificSites, ALLCOUNTED){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSGENCHR_getSpecificSites, EMPTY){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSGENCHR_getSpecificSites, THROWUGENE){
       ASSERT_TRUE(false);
 }
/**<  removeSpecificSites
 */
TEST(uBasicNGSGENCHR_removeSpecificSites, NONEREMOVED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_removeSpecificSites, SOMEREMOVED){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_removeSpecificSites, ALLREMOVED){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSGENCHR_removeSpecificSites, EMPTY){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSGENCHR_removeSpecificSites, THROWUGENE){
       ASSERT_TRUE(false);
 }

/**< applyOnAllSites */
TEST(uBasicNGSGENCHR_applyOnAllSites, EMPTY){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_applyOnAllSites, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_applyOnAllSites, EXCEPTION){
       ASSERT_TRUE(false);
 }
/**< applyOnAllSites */
TEST(uBasicNGSGENCHR_applyOnAllSitesConst, EMPTY){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_applyOnAllSitesConst, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_applyOnAllSitesConst, EXCEPTION){
       ASSERT_TRUE(false);
 }
/**< accumulateSitesInfo */
TEST(uBasicNGSGENCHR_accumulateSitesInfos, EMPTY){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_accumulateSitesInfo, NORMAL){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSGENCHR_accumulateSitesInfo, EXCEPTION){
       ASSERT_TRUE(false);
 }
/**< MinMaxSites */
TEST(uBasicNGSGENCHR_minAndMaxSites, NORMAL){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSGENCHR_minAndMaxSites, EMPTY){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSGENCHR_minAndMaxSites, EXCEPTION){
       ASSERT_TRUE(false);
 }

/**<countSitesWithProperty */
TEST(uBasicNGSGENCHR_countSitesWithProperty, NORMAL){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSGENCHR_countSitesWithProperty, EMPTY){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSGENCHR_countSitesWithProperty, EXCEPTION){
       ASSERT_TRUE(false);
 }
/**< applyOnAllChroms */


 TEST(uBasicNGSGENCHR_applyOnAllChroms, VALID){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSGENCHR_applyOnAllChroms, EMPTY){
       ASSERT_TRUE(false);
 }
