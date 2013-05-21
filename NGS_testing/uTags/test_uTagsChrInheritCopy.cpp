
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

class StandardChromsTags
{
public:
	StandardChromsTags()
	{
		/**< m_uTagsChroms[0]: Empty chrom without a name (""). */
       oneChr.setChr("chr1");

       oneChr.addData(uTags("chr1",100,200));

       manyChr.setChr("chr1");
       manyChr.addData(uTags("chr1",100,200));
       manyChr.addData(uTags("chr1",230,300));
       manyChr.addData(uTags("chr1",120,250));

       emptyChr.setChr("chr1");
	}
	uTagsChrom oneChr;
    uTagsChrom manyChr;
    uTagsChrom emptyChr;
};



























TEST(uTagsGENCHR_applyAndGetChrom, NORMAL){

       StandardChromsTags testChroms;
       auto functOp = [](uTags & item){   item.extendSite(20);
       };
       uTagsChrom results=testChroms.manyChr.applyAndGetChrom(functOp);

       EXPECT_TRUE(results.getSite(0).isEqual(uTags("chr1",80,220)));
       EXPECT_TRUE(results.getSite(1).isEqual(uTags("chr1",210,320)));
       EXPECT_TRUE(results.getSite(2).isEqual(uTags("chr1",100,270)));
 }
TEST(uTagsGENCHR_applyAndGetChrom, SIDEEFFECT){
       StandardChromsTags testChroms;
       long int counter=0;
       auto functOp = [&](uTags & item){   item.extendSite(20);
       counter+=20;
       };
       uTagsChrom results=testChroms.manyChr.applyAndGetChrom(functOp);
       EXPECT_TRUE(results.getSite(0).isEqual(uTags("chr1",80,220)));
       EXPECT_TRUE(results.getSite(1).isEqual(uTags("chr1",210,320)));
       EXPECT_TRUE(results.getSite(2).isEqual(uTags("chr1",100,270)));
       EXPECT_EQ(results.count(),3);
       EXPECT_EQ(counter,60);
 }
TEST(uTagsGENCHR_applyAndGetChrom, EMPTY){
       StandardChromsTags testChroms;
       auto functOp = [](uTags & item){   item.extendSite(20);
       };
       uTagsChrom results=testChroms.emptyChr.applyAndGetChrom(functOp);
       EXPECT_EQ(results.count(),0);
 }

/**< Apply on AllVecData */
TEST(uTagsGENCHR_applyAndGetVecData, NORMAL){
        StandardChromsTags testChroms;
       auto functOp = [&](uTags & item){   item.extendSite(20);
       };
       auto results=testChroms.manyChr.applyAndGetVecData(functOp);
       EXPECT_TRUE(results.at(0).isEqual(uTags("chr1",80,220)));
       EXPECT_TRUE(results.at(1).isEqual(uTags("chr1",210,320)));
       EXPECT_TRUE(results.at(2).isEqual(uTags("chr1",100,270)));

 }
TEST(uTagsGENCHR_applyAndGetVecData, SIDEEFFECT){
       StandardChromsTags testChroms;
       long int counter=0;
       auto functOp = [&](uTags & item){   item.extendSite(20);
       counter+=20;
       };
       auto results=testChroms.manyChr.applyAndGetVecData(functOp);
       EXPECT_TRUE(results.at(0).isEqual(uTags("chr1",80,220)));
       EXPECT_TRUE(results.at(1).isEqual(uTags("chr1",210,320)));
       EXPECT_TRUE(results.at(2).isEqual(uTags("chr1",100,270)));
       EXPECT_EQ((int)results.size(),3);
       EXPECT_EQ(counter,60);
 }
TEST(uTagsGENCHR_applyAndGetVecData, EMPTY){
       StandardChromsTags testChroms;
       auto functOp = [](uTags & item){   item.extendSite(20);
       };
       auto results=testChroms.emptyChr.applyAndGetVecData(functOp);
       EXPECT_EQ((int)results.size(),0);
 }

/**<  computeOnAllSites*/
TEST(uTagsGENCHR_computeOnAllSites, NORMAL){
       StandardChromsTags testChroms;
       auto functOp = [&](uTags item)->int{  return item.getLength();};
       auto results=testChroms.manyChr.computeOnAllSites(functOp);
       EXPECT_EQ(results.at(0),testChroms.manyChr.getSite(0).getLength());
       EXPECT_EQ(results.at(1),testChroms.manyChr.getSite(1).getLength());
       EXPECT_EQ(results.at(2),testChroms.manyChr.getSite(2).getLength());
       EXPECT_EQ((int)results.size(),testChroms.manyChr.count());

       EXPECT_EQ(std::accumulate(results.begin(), results.end(), 0), (int)testChroms.manyChr.sumSiteSize());
 }
TEST(uTagsGENCHR_computeOnAllSites, EMPTY){
       StandardChromsTags testChroms;
       auto functOp = [&](uTags item)->int{  return item.getLength();};
       auto results=testChroms.emptyChr.computeOnAllSites(functOp);
       EXPECT_EQ((int)results.size(),testChroms.emptyChr.count());
 }

/**<  getSpecificSites*/
TEST(uTagsGENCHR_getSpecificSites, NONECOUNTED){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (item.getLength()>2000);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),0);
 }
TEST(uTagsGENCHR_getSpecificSites, SOMECOUNTED){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (item.getLength()>99);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),2);
 }
TEST(uTagsGENCHR_getSpecificSites, ALLCOUNTED){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (item.getLength()>5);
       };
       auto results=testChroms.manyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),3);
 }
 TEST(uTagsGENCHR_getSpecificSites, EMPTY){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (true);
       };
       auto results=testChroms.emptyChr.getSpecificSites(functOp);
       EXPECT_EQ((int)results.size(),0);
 }
/**<  removeSpecificSites
 */
TEST(uTagsGENCHR_removeSpecificSites, NONEREMOVED){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (item.getLength()>2000); };

       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),3);
 }
TEST(uTagsGENCHR_removeSpecificSites, SOMEREMOVED){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (item.getLength()>99);
       };
       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),1);
 }
TEST(uTagsGENCHR_removeSpecificSites, ALLREMOVED){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (item.getLength()>5);
       };
       testChroms.manyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.manyChr.count(),0);
 }
 TEST(uTagsGENCHR_removeSpecificSites, EMPTY){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (false);
       };
       testChroms.emptyChr.removeSpecificSites(functOp);
       EXPECT_EQ(testChroms.emptyChr.count(),0);
 }

/**< applyOnAllSites */
TEST(uTagsGENCHR_applyOnAllSites, EMPTY){

       StandardChromsTags testChroms;
       auto functOp = [](uTags & item){   item.extendSite(20);
       };
       EXPECT_NO_THROW(testChroms.manyChr.applyOnAllSites(functOp));
 }
TEST(uTagsGENCHR_applyOnAllSites, NORMAL){
        StandardChromsTags testChroms;
       auto functOp = [](uTags & item){   item.extendSite(20);
       };
       testChroms.manyChr.applyOnAllSites(functOp);
       EXPECT_TRUE(testChroms.manyChr.getSite(0).isEqual(uTags("chr1",80,220)));
       EXPECT_TRUE(testChroms.manyChr.getSite(1).isEqual(uTags("chr1",210,320)));
       EXPECT_TRUE(testChroms.manyChr.getSite(2).isEqual(uTags("chr1",100,270)));
 }
TEST(uTagsGENCHR_applyOnAllSites, EXCEPTION){

    StandardChromsTags testChroms;
       auto functOp = [](uTags & item){   throw param_throw();
       };

      EXPECT_THROW(testChroms.manyChr.applyOnAllSites(functOp),param_throw);

 }
/**< applyOnAllSites */
TEST(uTagsGENCHR_applyOnAllSitesConst, EMPTY){
       const StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){
       };
       EXPECT_NO_THROW(testChroms.manyChr.applyOnAllSites(functOp));
 }
TEST(uTagsGENCHR_applyOnAllSitesConst, NORMAL){
       const StandardChromsTags testChroms;
       int siteCount=0;
       auto functOp = [&](const uTags & item){siteCount+=item.getLength();
       };
       testChroms.manyChr.applyOnAllSites(functOp);
       EXPECT_EQ(siteCount,  testChroms.manyChr.sumSiteSize());
 }
TEST(uTagsGENCHR_applyOnAllSitesConst, EXCEPTION){
 StandardChromsTags testChroms;
       auto functOp = [](uTags & item){   throw param_throw();
       };
      EXPECT_THROW(testChroms.manyChr.applyOnAllSites(functOp),param_throw);

 }
/**< accumulateSitesInfo */
TEST(uTagsGENCHR_accumulateSitesInfos, EMPTY){
       const StandardChromsTags testChroms;
       int siteCount=0;
       auto functOp = [&](int siteCounts,const uTags & item){ return (siteCounts+=item.getLength());
       };
       EXPECT_EQ(testChroms.emptyChr.accumulateSitesInfo(functOp,siteCount),  testChroms.emptyChr.sumSiteSize());
 }
TEST(uTagsGENCHR_accumulateSitesInfo, NORMAL){
       const StandardChromsTags testChroms;
       int siteCount=0;
       auto functOp = [&](int siteCounts,const uTags & item){return siteCounts+=item.getLength();
       };
       EXPECT_EQ(testChroms.manyChr.accumulateSitesInfo(functOp,siteCount),  testChroms.manyChr.sumSiteSize());
 }
TEST(uTagsGENCHR_accumulateSitesInfo, EXCEPTION){

       const StandardChromsTags testChroms;
       auto functOp = [](int siteCounts,const uTags & item){   throw param_throw();
       return 0;
       };
       int siteCount =0;
       EXPECT_THROW(testChroms.manyChr.accumulateSitesInfo(functOp,siteCount),param_throw);
 }


/**<countSitesWithProperty */
TEST(uTagsGENCHR_countSitesWithProperty, SOMECOUNTED){
      StandardChromsTags testChroms;
       auto functOp = [](const uTags & item) mutable {  return (item.getLength()>99);
       StandardChromsTags testChroms;
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),2);
 }

 TEST(uTagsGENCHR_countSitesWithProperty, ALLCOUNTED){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (item.getLength()>2);
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),3);
 }

 TEST(uTagsGENCHR_countSitesWithProperty, NONECOUNTED){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return (item.getLength()>200000);
       };
       EXPECT_EQ(  testChroms.manyChr.countSitesWithProperty(functOp),0);
 }
  TEST(uTagsGENCHR_countSitesWithProperty, EMPTY){
       StandardChromsTags testChroms;
       auto functOp = [](const uTags & item){  return true;
       };
       EXPECT_EQ(  testChroms.emptyChr.countSitesWithProperty(functOp),0);
 }

/**< minSite(comp ) */
TEST(uTagsGENCHR_minSite, NORMAL)
{
       StandardChromsTags testChroms;
       auto minItr=  testChroms.manyChr.minSite(uTagsChrom::comparePos);
       EXPECT_EQ(minItr->getStart(),uTags("chr1",100,200).getStart() );
}
TEST(uTagsGENCHR_minSite, EMPTY)
{
       StandardChromsTags testChroms;
       auto minItr=  testChroms.emptyChr.minSite(uTagsChrom::comparePos);
       EXPECT_EQ(minItr, testChroms.emptyChr.end() );
}
TEST(uTagsGENCHR_minSite, CUSTOM)
{
       StandardChromsTags testChroms;
       auto minItr=  testChroms.manyChr.minSite(uTagsChrom::compareLength);
       EXPECT_EQ(minItr->getStart(), uTags("chr1",230,300).getStart() );
}
/**< maxSite(comp ) */

TEST(uTagsGENCHR_maxSite, NORMAL)
{
       StandardChromsTags testChroms;
       auto maxItr=  testChroms.manyChr.maxSite(uTagsChrom::comparePos);
       EXPECT_EQ(maxItr->getStart(),uTags("chr1",230,300).getStart() );
}
TEST(uTagsGENCHR_maxSite, EMPTY)
{
       StandardChromsTags testChroms;
       auto minItr=  testChroms.emptyChr.maxSite(uTagsChrom::comparePos);
       EXPECT_EQ(minItr, testChroms.emptyChr.end() );
}
TEST(uTagsGENCHR_maxSite, CUSTOM)
{
       StandardChromsTags testChroms;
       auto minItr=  testChroms.manyChr.maxSite(uTagsChrom::compareLength);
       EXPECT_EQ(minItr->getStart(), uTags("chr1",120,250).getStart() );
}
/**< MinMaxSites */
TEST(uTagsGENCHR_minAndMaxSites, VALIDCOMPILE){
       StandardChromsTags testChroms;
        testChroms.manyChr.minAndMaxSites(uTagsChrom::compareLength);
 }
