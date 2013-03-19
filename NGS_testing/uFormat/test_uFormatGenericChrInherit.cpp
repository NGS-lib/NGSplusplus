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

class StandardChroms
{
public:
	StandardChroms()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
       oneChr.setChr("chr1");

       oneChr.addData(uBasicNGS("chr1",100,200));

       manyChr.setChr("chr1");
       manyChr.addData(uBasicNGS("chr1",100,200));
       manyChr.addData(uBasicNGS("chr1",230,300));
       manyChr.addData(uBasicNGS("chr1",120,250));

       emptyChr.setChr("chr1");
	}
	uBasicNGSChrom oneChr;
    uBasicNGSChrom manyChr;
    uBasicNGSChrom emptyChr;
};


TEST(uBasicNGSGENCHR_applyAndGetChrom, NORMAL){

       StandardChroms testChroms;
       auto functOp = [](uBasicNGS & item){   item.extendSite(20);
       };
       uBasicNGSChrom results=testChroms.manyChr.applyAndGetChrom(functOp);
       EXPECT_TRUE(results.getSite(0).isEqual(uBasicNGS("chr1",80,220)));
       EXPECT_TRUE(results.getSite(1).isEqual(uBasicNGS("chr1",210,320)));
       EXPECT_TRUE(results.getSite(2).isEqual(uBasicNGS("chr1",100,270)));
 }
TEST(uBasicNGSGENCHR_applyAndGetChrom, SIDEEFFECT){
       StandardChroms testChroms;
       long int counter=0;
       auto functOp = [&](uBasicNGS & item){   item.extendSite(20);
       counter+=20;
       };
       uBasicNGSChrom results=testChroms.manyChr.applyAndGetChrom(functOp);
       EXPECT_TRUE(results.getSite(0).isEqual(uBasicNGS("chr1",80,220)));
       EXPECT_TRUE(results.getSite(1).isEqual(uBasicNGS("chr1",210,320)));
       EXPECT_TRUE(results.getSite(2).isEqual(uBasicNGS("chr1",100,270)));
       EXPECT_EQ(results.count(),3);
       EXPECT_EQ(counter,60);
 }
TEST(uBasicNGSGENCHR_applyAndGetChrom, EMPTY){
       StandardChroms testChroms;
       auto functOp = [](uBasicNGS & item){   item.extendSite(20);
       };
       uBasicNGSChrom results=testChroms.emptyChr.applyAndGetChrom(functOp);
       EXPECT_EQ(results.count(),0);
 }

/**< Apply on AllVecData */
TEST(uBasicNGSGENCHR_applyAndGetVecData, NORMAL){
        StandardChroms testChroms;
       auto functOp = [&](uBasicNGS & item){   item.extendSite(20);
       };
       auto results=testChroms.manyChr.applyAndGetVecData(functOp);
       EXPECT_TRUE(results.at(0).isEqual(uBasicNGS("chr1",80,220)));
       EXPECT_TRUE(results.at(1).isEqual(uBasicNGS("chr1",210,320)));
       EXPECT_TRUE(results.at(2).isEqual(uBasicNGS("chr1",100,270)));

 }
TEST(uBasicNGSGENCHR_applyAndGetVecData, SIDEEFFECT){
       StandardChroms testChroms;
       long int counter=0;
       auto functOp = [&](uBasicNGS & item){   item.extendSite(20);
       counter+=20;
       };
       auto results=testChroms.manyChr.applyAndGetVecData(functOp);
       EXPECT_TRUE(results.at(0).isEqual(uBasicNGS("chr1",80,220)));
       EXPECT_TRUE(results.at(1).isEqual(uBasicNGS("chr1",210,320)));
       EXPECT_TRUE(results.at(2).isEqual(uBasicNGS("chr1",100,270)));
       EXPECT_EQ(results.size(),3);
       EXPECT_EQ(counter,60);
 }
TEST(uBasicNGSGENCHR_applyAndGetVecData, EMPTY){
       StandardChroms testChroms;
       auto functOp = [](uBasicNGS & item){   item.extendSite(20);
       };
       auto results=testChroms.emptyChr.applyAndGetVecData(functOp);
       EXPECT_EQ(results.size(),0);
 }

/**<  computeOnAllSites*/
TEST(uBasicNGSGENCHR_computeOnAllSites, NORMAL){
       StandardChroms testChroms;
       auto functOp = [&](uBasicNGS item)->int{  return item.getLenght();};
       auto results=testChroms.manyChr.computeOnAllSites(functOp);
       EXPECT_EQ(results.at(0),testChroms.manyChr.getSite(0).getLenght());
       EXPECT_EQ(results.at(1),testChroms.manyChr.getSite(1).getLenght());
       EXPECT_EQ(results.at(2),testChroms.manyChr.getSite(2).getLenght());
       EXPECT_EQ(results.size(),testChroms.manyChr.count());

       EXPECT_EQ(std::accumulate(results.begin(), results.end(), 0), testChroms.manyChr.sumSiteSize());
 }
TEST(uBasicNGSGENCHR_computeOnAllSites, EMPTY){
       StandardChroms testChroms;
       auto functOp = [&](uBasicNGS item)->int{  return item.getLenght();};
       auto results=testChroms.emptyChr.computeOnAllSites(functOp);
       EXPECT_EQ(results.size(),testChroms.emptyChr.count());
 }

/**<  getSpecificSites*/
TEST(uBasicNGSGENCHR_getSpecificSites, NONECOUNTED){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>2000);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ(results.size(),0);
 }
TEST(uBasicNGSGENCHR_getSpecificSites, SOMECOUNTED){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>99);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ(results.size(),2);
 }
TEST(uBasicNGSGENCHR_getSpecificSites, ALLCOUNTED){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>5);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ(results.size(),3);
 }
 TEST(uBasicNGSGENCHR_getSpecificSites, EMPTY){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (true);
       };
       auto results=testChroms.emptyChr.getSpecificSites(functOp);
       EXPECT_EQ(results.size(),0);
 }
/**<  removeSpecificSites
 */
TEST(uBasicNGSGENCHR_removeSpecificSites, NONEREMOVED){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>2000); };

       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),3);
 }
TEST(uBasicNGSGENCHR_removeSpecificSites, SOMEREMOVED){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>99);
       };
       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),1);
 }
TEST(uBasicNGSGENCHR_removeSpecificSites, ALLREMOVED){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>5);
       };
       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),0);
 }
 TEST(uBasicNGSGENCHR_removeSpecificSites, EMPTY){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (false);
       };
       testChroms.emptyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.emptyChr.count(),0);
 }

/**< applyOnAllSites */
TEST(uBasicNGSGENCHR_applyOnAllSites, EMPTY){

       StandardChroms testChroms;
       auto functOp = [](uBasicNGS & item){   item.extendSite(20);
       };
       EXPECT_NO_THROW(testChroms.manyChr.applyOnAllSites(functOp));
 }
TEST(uBasicNGSGENCHR_applyOnAllSites, NORMAL){
        StandardChroms testChroms;
       auto functOp = [](uBasicNGS & item){   item.extendSite(20);
       };
       testChroms.manyChr.applyOnAllSites(functOp);
       EXPECT_TRUE(testChroms.manyChr.getSite(0).isEqual(uBasicNGS("chr1",80,220)));
       EXPECT_TRUE(testChroms.manyChr.getSite(1).isEqual(uBasicNGS("chr1",210,320)));
       EXPECT_TRUE(testChroms.manyChr.getSite(2).isEqual(uBasicNGS("chr1",100,270)));
 }
TEST(uBasicNGSGENCHR_applyOnAllSites, EXCEPTION){

    StandardChroms testChroms;
       auto functOp = [](uBasicNGS & item){   throw param_throw();
       };

      EXPECT_THROW(testChroms.manyChr.applyOnAllSites(functOp),param_throw);

 }
/**< applyOnAllSites */
TEST(uBasicNGSGENCHR_applyOnAllSitesConst, EMPTY){
       const StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){
       };
       EXPECT_NO_THROW(testChroms.manyChr.applyOnAllSites(functOp));
 }
TEST(uBasicNGSGENCHR_applyOnAllSitesConst, NORMAL){
       const StandardChroms testChroms;
       int siteCount=0;
       auto functOp = [&](const uBasicNGS & item){siteCount+=item.getLenght();
       };
       testChroms.manyChr.applyOnAllSites(functOp);
       EXPECT_EQ(siteCount,  testChroms.manyChr.sumSiteSize());
 }
TEST(uBasicNGSGENCHR_applyOnAllSitesConst, EXCEPTION){
 StandardChroms testChroms;
       auto functOp = [](uBasicNGS & item){   throw param_throw();
       };
      EXPECT_THROW(testChroms.manyChr.applyOnAllSites(functOp),param_throw);

 }
/**< accumulateSitesInfo */
TEST(uBasicNGSGENCHR_accumulateSitesInfos, EMPTY){
       const StandardChroms testChroms;
       int siteCount=0;
       auto functOp = [&](int siteCounts,const uBasicNGS & item){ return (siteCounts+=item.getLenght());
       };
       EXPECT_EQ(testChroms.emptyChr.accumulateSitesInfo(functOp,siteCount),  testChroms.emptyChr.sumSiteSize());
 }
TEST(uBasicNGSGENCHR_accumulateSitesInfo, NORMAL){
       const StandardChroms testChroms;
       int siteCount=0;
       auto functOp = [&](int siteCounts,const uBasicNGS & item){return siteCounts+=item.getLenght();
       };
       EXPECT_EQ(testChroms.manyChr.accumulateSitesInfo(functOp,siteCount),  testChroms.manyChr.sumSiteSize());
 }
TEST(uBasicNGSGENCHR_accumulateSitesInfo, EXCEPTION){

       const StandardChroms testChroms;
       auto functOp = [](int siteCounts,const uBasicNGS & item){   throw param_throw();
       return 0;
       };
       int siteCount =0;
       EXPECT_THROW(testChroms.manyChr.accumulateSitesInfo(functOp,siteCount),param_throw);
 }


/**<countSitesWithProperty */
TEST(uBasicNGSGENCHR_countSitesWithProperty, SOMECOUNTED){
      StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item) mutable {  return (item.getLenght()>99);
       StandardChroms testChroms;
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),2);
 }

 TEST(uBasicNGSGENCHR_countSitesWithProperty, ALLCOUNTED){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>2);
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),3);
 }

 TEST(uBasicNGSGENCHR_countSitesWithProperty, NONECOUNTED){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>200000);
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),0);
 }
  TEST(uBasicNGSGENCHR_countSitesWithProperty, EMPTY){
       StandardChroms testChroms;
       auto functOp = [](const uBasicNGS & item){  return true;
       };
       EXPECT_EQ(  testChroms.emptyChr.countSitesWithProperty(functOp),0);
 }

/**< minSite(comp ) */
TEST(uBasicNGSGENCHR_minSite, NORMAL)
{
       StandardChroms testChroms;
       auto minItr=  testChroms.manyChr.minSite(uBasicNGSChrom::comparePos);
       EXPECT_EQ(minItr->getStart(),uBasicNGS("chr1",100,200).getStart() );
}
TEST(uBasicNGSGENCHR_minSite, EMPTY)
{
       StandardChroms testChroms;
       auto minItr=  testChroms.emptyChr.minSite(uBasicNGSChrom::comparePos);
       EXPECT_EQ(minItr, testChroms.emptyChr.end() );
}
TEST(uBasicNGSGENCHR_minSite, CUSTOM)
{
       StandardChroms testChroms;
       auto minItr=  testChroms.manyChr.minSite(uBasicNGSChrom::compareLenght);
       EXPECT_EQ(minItr->getStart(), uBasicNGS("chr1",230,300).getStart() );
}
/**< maxSite(comp ) */

TEST(uBasicNGSGENCHR_maxSite, NORMAL)
{
       StandardChroms testChroms;
       auto maxItr=  testChroms.manyChr.maxSite(uBasicNGSChrom::comparePos);
       EXPECT_EQ(maxItr->getStart(),uBasicNGS("chr1",230,300).getStart() );
}
TEST(uBasicNGSGENCHR_maxSite, EMPTY)
{

       StandardChroms testChroms;
       auto minItr=  testChroms.emptyChr.maxSite(uBasicNGSChrom::comparePos);
       EXPECT_EQ(minItr, testChroms.emptyChr.end() );
}
TEST(uBasicNGSGENCHR_maxSite, CUSTOM)
{
       StandardChroms testChroms;
       auto minItr=  testChroms.manyChr.maxSite(uBasicNGSChrom::compareLenght);
       EXPECT_EQ(minItr->getStart(), uBasicNGS("chr1",120,250).getStart() );
}
/**< MinMaxSites */
// TODO: more tests
TEST(uBasicNGSGENCHR_minAndMaxSites, VALIDCOMPILE){
       StandardChroms testChroms;
        testChroms.manyChr.minAndMaxSites(uBasicNGSChrom::compareLenght);
 }

