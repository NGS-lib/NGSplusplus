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

/**< FIXTURES: Begin */
class validChroms
{
public:
	validChroms()
	{
		/**< m_uBasicNGSChroms[0]: Empty chrom without a name (""). */
		uBasicNGSChrom chrom_0;
		chrom_0.setChr("");
		m_uBasicNGSChroms.push_back(chrom_0);

		/**< m_uBasicNGSChroms[1]: Chrom with 1 element without a name ("") */
		uBasicNGSChrom chrom_1;
		chrom_1.setChr("");
		chrom_1.addData(uBasicNGS("",100,200));
		m_uBasicNGSChroms.push_back(chrom_1);
		
		/**< m_uBasicNGSChroms[2]: Chrom with 3 elements without a name ("") */
		uBasicNGSChrom chrom_2;
		chrom_2.setChr("");
		chrom_2.addData(uBasicNGS("",100,200));
		chrom_2.addData(uBasicNGS("",200,300));
		chrom_2.addData(uBasicNGS("",300,400));
		m_uBasicNGSChroms.push_back(chrom_2);

		/**< m_uBasicNGSChroms[3]: Empty chrom with a name ("chr3") */
		uBasicNGSChrom chrom_3;
		chrom_3.setChr("chr3");
		m_uBasicNGSChroms.push_back(chrom_3);

		/**< m_uBasicNGSChroms[4]: Chrom with 1 element with a name ("chr4") */
		uBasicNGSChrom chrom_4;
		chrom_4.setChr("chr4");
		chrom_4.addData(uBasicNGS("chr4",100,200));
		m_uBasicNGSChroms.push_back(chrom_4);

		/**< m_uBasicNGSChroms[5]: Chrom with 3 elements with a name and no overlapping element ("chr5") */
		uBasicNGSChrom chrom_5;
		chrom_5.setChr("chr5");
		chrom_5.addData(uBasicNGS("chr5",100,200));
		chrom_5.addData(uBasicNGS("chr5",300,400));
		chrom_5.addData(uBasicNGS("chr5",500,600));
		m_uBasicNGSChroms.push_back(chrom_5);

		/**< m_uBasicNGSChroms[6]: Chrom with 3 elements with a name and 2 overlapping elements ("chr6") */
		uBasicNGSChrom chrom_6;
		chrom_6.setChr("chr6");
		chrom_6.addData(uBasicNGS("chr6",100,250));
		chrom_6.addData(uBasicNGS("chr6",200,300));
		chrom_6.addData(uBasicNGS("chr6",500,600));
		m_uBasicNGSChroms.push_back(chrom_6);

		/**< m_uBasicNGSChroms[7]: Chrom with 3 elements with a name and 3 overlapping elements ("chr7") */
		uBasicNGSChrom chrom_7;
		chrom_7.setChr("chr7");
		chrom_7.addData(uBasicNGS("chr7",100,250));
		chrom_7.addData(uBasicNGS("chr7",200,350));
		chrom_7.addData(uBasicNGS("chr7",300,400));
		m_uBasicNGSChroms.push_back(chrom_7);
	}
	vector<uBasicNGSChrom> m_uBasicNGSChroms;
};

//class test_uFormatExpInherit_multipleChrom : public testing::Test {
//public:
//	test_uFormatExpInherit_multipleChrom() 
//	{
//
//	}
// protected:
///**< As always, this is inclusive so 100-199 is of size 100 */
//  virtual void SetUp() {
//
//    uChromTestOverlap.setChr("chr1");
//    uChromTestOverlap.addData(uBasicNGS ("chr1", 300, CHROMDIVIDESIZE));
//    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 199));
//    uChromTestOverlap.addData(uBasicNGS("chr1", 100, 299));
//  }
//uBasicNGSChrom uChromTestOverlap;
//};


/**< FIXTURES: End */

/*
 * Test for the funciton:
 *		_CHROM_ getChrom(const std::string & chrom) const;
 *	Valid cases:
 *		NONAMECHROM
 *		VALIDCHROM
 *	Invalid cases:
 *		NOCHROMTHROWEXC
 */

TEST(uBasicNGSEXP_GetChrom, NONAMECHROM)
{
 	uBasicNGSExperiment anExp;
	validChroms chroms;
	anExp.addData(chroms.m_uBasicNGSChroms[1]); // Chrom with no name and 1 element
	uBasicNGSChrom aChrom = anExp.getChrom("");
	EXPECT_EQ(aChrom.count(), 1);
	EXPECT_EQ(aChrom.getSite(0).getChr(), "");
	EXPECT_EQ(aChrom.getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom.getSite(0).getEnd(), 200);
}

TEST(uBasicNGSEXP_GetChrom, VALIDCHROM)
{
 	uBasicNGSExperiment anExp;
	validChroms chroms;
	anExp.addData(chroms.m_uBasicNGSChroms[4]); // Chrom with a name and 1 element
	uBasicNGSChrom aChrom = anExp.getChrom("chr4");
	EXPECT_EQ(aChrom.count(), 1);
	EXPECT_EQ(aChrom.getSite(0).getChr(), "chr4");
	EXPECT_EQ(aChrom.getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom.getSite(0).getEnd(), 200);
}

TEST(uBasicNGSEXP_GetChrom, NOCHROMTHROWEXC)
{
	uBasicNGSExperiment anExp;
	EXPECT_THROW(anExp.getChrom(""), ugene_operation_throw);
}

/**<  */
TEST(uBasicNGSEXP_getChromP, NONAMECHROM){
    //Chromosome pointer no name and works
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getChromP, NOCHROMTHROWEXC){
    //Chromosome pointer not exist and throw
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getChromP, VALIDCHROM){
    // Chromosome pointer, valid and works
       ASSERT_TRUE(false);
 }
/**<  */
 TEST(uBasicNGSEXP_getConstChromP, NONAMECHROM){
    //Chromosome pointer no name and works
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getConstChromP, NOCHROMTHROWEXC){
    //Chromosome pointer not exist and throw
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getConstChromP, VALIDCHROM){
    // Chromosome pointer, valid and works
       ASSERT_TRUE(false);
 }
/**<  */

/**< GetSite */
TEST(uBasicNGSEXP_getSite, VALID){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getSite, OUTOFBOUND){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSite, BELOW0){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSite, VALIDITERRATOR){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSite, INVALIDITERATOR){
       ASSERT_TRUE(false);
 }
 /**< GetOverlapping */
 TEST(uBasicNGSEXP_getOverlapping, VALIDEXP){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, VALIDCHROM){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, VALIDPOS){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getOverlapping, EMPTYEXP){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, EMPTYCHROM){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISEXP){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISCHROM){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISANDPOS){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getOverlapping, POLYMORPHICHROM){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getOverlapping, POLYMORPHICEXP){
       ASSERT_TRUE(false);
 }

/**< set/get Chr Size */
TEST(uBasicNGSEXP_getChrSize, VALID){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getChrSize, NOTVALID){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getChrSize, EMPTY){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSEXP_setOverlapping, INVALIDVALUE){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_setOverlapping, INVALUECHR){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_setOverlapping, VALID){
       ASSERT_TRUE(false);
 }

/**< getSubsetCount */
 TEST(uBasicNGSEXP_getSubsetCount, POSITIONS){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getSubsetCount, ELEMENT){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSubsetCount, CHROMNOEXIST){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSubsetCount, NONAMECHROM){
       ASSERT_TRUE(false);
 }
