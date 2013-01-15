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

/**< pseudo-FIXTURES: Begin */
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

class validExperiments
{
public:
	validExperiments()
	{
		validChroms chroms;
		/**< NoName_Empty Chrom */
		uBasicNGSExperiment NoName_EmptyExp;
		NoName_EmptyExp.addData(chroms.m_uBasicNGSChroms[0]);
		m_uBasicNGSExp["NoName_Empty"] = NoName_EmptyExp;

		/**< NoName_1elem Chrom */
		uBasicNGSExperiment NoName_1elemExp;
		NoName_1elemExp.addData(chroms.m_uBasicNGSChroms[1]);
		m_uBasicNGSExp["NoName_1elem"] = NoName_1elemExp;

		/**< NoName_3elems Chrom */
		uBasicNGSExperiment NoName_3elemsExp;
		NoName_3elemsExp.addData(chroms.m_uBasicNGSChroms[2]);
		m_uBasicNGSExp["NoName_3elems"] = NoName_3elemsExp;

		/**< MultipleChroms */
		uBasicNGSExperiment MultipleChromsExp;
		for (size_t i = 2; i < chroms.m_uBasicNGSChroms.size(); i++)
		{
			MultipleChromsExp.addData(chroms.m_uBasicNGSChroms[i]); // Chrom with no name and 1 element
		}
		m_uBasicNGSExp["MultipleChroms"] = MultipleChromsExp;

		/**< Empty Exp */
		uBasicNGSExperiment Empty_Exp;
		m_uBasicNGSExp["Empty_Exp"] = Empty_Exp;

	}
	uBasicNGSExperiment* getExperiment(const std::string& name)
	{
		return &(m_uBasicNGSExp[name]);
	}

	map<std::string, uBasicNGSExperiment> m_uBasicNGSExp;
};
/**< pseudo-FIXTURES: End */

/*
 * Test for the function:
 *		_CHROM_ getChrom(const std::string & chrom) const;
 *	Valid cases:
 *		NONAMECHROM
 *		VALIDCHROM
 *	Invalid cases:
 *		NOCHROMTHROWEXC
 */

TEST(uBasicNGSEXP_GetChrom, NONAMECHROM)
{
	validExperiments myExperiments;
	uBasicNGSChrom aChrom = myExperiments.getExperiment("NoName_1elem")->getChrom("");
	EXPECT_EQ(aChrom.count(), 1);
	EXPECT_EQ(aChrom.getSite(0).getChr(), "");
	EXPECT_EQ(aChrom.getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom.getSite(0).getEnd(), 200);
}

TEST(uBasicNGSEXP_GetChrom, VALIDCHROM)
{
	validExperiments myExperiments;
	uBasicNGSChrom aChrom = myExperiments.getExperiment("MultipleChroms")->getChrom("chr4");
	EXPECT_EQ(aChrom.count(), 1);
	EXPECT_EQ(aChrom.getSite(0).getChr(), "chr4");
	EXPECT_EQ(aChrom.getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom.getSite(0).getEnd(), 200);
}

TEST(uBasicNGSEXP_GetChrom, NOCHROMTHROWEXC)
{
	uBasicNGSExperiment myExperiments;
	EXPECT_THROW(myExperiments.getChrom(""), ugene_operation_throw);
}

/*
 * Test for the function:
 *		_CHROM_* getpChrom(const std::string & chrom);
 *	Valid cases:
 *		NONAMECHROM
 *		VALIDCHROM
 *	Invalid cases:
 *		NOCHROMTHROWEXC
 */
 
TEST(uBasicNGSEXP_getpChrom, NONAMECHROM)
{
	validExperiments myExperiments;
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("NoName_1elem")->getpChrom("");
	EXPECT_EQ(aChrom->count(), 1);
	EXPECT_EQ(aChrom->getSite(0).getChr(), "");
	EXPECT_EQ(aChrom->getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom->getSite(0).getEnd(), 200);
//	EXPECT_NO_THROW(aChrom->addData(uBasicNGS("chr4", 1000, 2000))); // TODO: Why does this fail?
}

TEST(uBasicNGSEXP_getpChrom, VALIDCHROM)
{
	validExperiments myExperiments;
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4");
	EXPECT_EQ(aChrom->count(), 1);
	EXPECT_EQ(aChrom->getSite(0).getChr(), "chr4");
	EXPECT_EQ(aChrom->getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom->getSite(0).getEnd(), 200);
//	EXPECT_NO_THROW(aChrom->addData(uBasicNGS("chr4", 1000, 2000))); // TODO: Why does this fail?
}

TEST(uBasicNGSEXP_getpChrom, NOCHROMTHROWEXC)
{
	uBasicNGSExperiment anExp;
	EXPECT_THROW(anExp.getpChrom(""), ugene_operation_throw);
}

/*
 * Test for the function:
 *		const _CHROM_* getpChrom(const std::string & chrom) const;
 *	Valid cases:
 *		NONAMECHROM
 *		VALIDCHROM
 *	Invalid cases:
 *		NOCHROMTHROWEXC
 */

TEST(uBasicNGSEXP_getpChromConst, NONAMECHROM)
{
	validExperiments myExperiments;
	const uBasicNGSChrom* aChrom = myExperiments.getExperiment("NoName_1elem")->getpChrom("");
	EXPECT_EQ((int)(aChrom->avgSiteSize()), 101);
	EXPECT_EQ(aChrom->count(), 1);
	EXPECT_EQ(aChrom->getSite(0).getChr(), "");
	EXPECT_EQ(aChrom->getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom->getSite(0).getEnd(), 200);
	// TODO: I don't see a way to test if it is impossible to call a non-const function since it won't compile.
//	EXPECT_THROW(aChrom->addData(uBasicNGS("chr4", 1000, 2000)), ugene_exception_base);
}

TEST(uBasicNGSEXP_getpChromConst, VALIDCHROM)
{
	validExperiments myExperiments;
	const uBasicNGSChrom* aChrom = myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4");
	EXPECT_EQ((int)(aChrom->avgSiteSize()), 101);
	EXPECT_EQ(aChrom->count(), 1);
	EXPECT_EQ(aChrom->getSite(0).getChr(), "chr4");
	EXPECT_EQ(aChrom->getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom->getSite(0).getEnd(), 200);
	// TODO: I don't see a way to test if it is impossible to call a non-const function since it won't compile.
//	EXPECT_THROW(aChrom->addData(uBasicNGS("chr4", 1000, 2000)), ugene_exception_base);
}

TEST(uBasicNGSEXP_getpChromConst, NOCHROMTHROWEXC)
{
	uBasicNGSExperiment anExp;
	try
	{
		const uBasicNGSChrom* aChrom = anExp.getpChrom("chr4");
		EXPECT_TRUE(false);
		aChrom->avgSiteSize(); /**< To avoid warning: unused variable */
	}
	catch (ugene_exception_base& e)
	{
		EXPECT_TRUE(true);
	}
}

/*
 * Test for the function:
 *		_BASE_ getSite(std::string chr, int position)const;
 *		_BASE_ getSite(typename std::vector<_BASE_>::const_iterator posItr)const;
 *	Valid cases:
 *		VALID
 *		VALIDITERRATOR
 *	Invalid cases:
 *		OUTOFBOUND
 *		BELOW0
 */

TEST(uBasicNGSEXP_getSite, VALID){
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSite("chr4", 0).getChr(), "chr4");
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSite("chr4", 0).getStart(), 100);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSite("chr4", 0).getEnd(), 200);
}

TEST(uBasicNGSEXP_getSite, VALIDITERRATOR){
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	auto it = myExperiments.getExperiment("MultipleChroms")->findNextSite("chr5", 0);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSite(it).getChr(), "chr5");
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSite(it).getStart(), 100); 
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSite(it).getEnd(),  200);
}

TEST(uBasicNGSEXP_getSite, OUTOFBOUND){
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->getSite("chr4", 1), std::exception);
}

TEST(uBasicNGSEXP_getSite, BELOW0){
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->getSite("chr4", -1), std::exception);
}


/*
 * Test for the function:
 *		_SELF_ getOverlapping(uGenericNGSExperiment<_SELFPAR_, _CHROMPAR_,_BASEPAR_> &compareExp, OverlapType type=OverlapType::OVERLAP_PARTIAL);
 *		template<class _SELFPAR_,typename _BASEPAR_>
 *			_SELF_ getOverlapping(uGenericNGSChrom<_SELFPAR_,_BASEPAR_> &compareExp, OverlapType type=OverlapType::OVERLAP_PARTIAL);
 *		_SELF_ getOverlapping(std::string chr, int start, int end, OverlapType type=OverlapType::OVERLAP_PARTIAL);
 *	Valid cases:
 *		VALIDEXP
 *		VALIDCHROM
 *		VALIDPOS
 *		EMPTYEXP
 *		EMPTYCHROM
 *		EMPTYTHISEXP
 *		EMPTYTHISCHROM
 *		EMPTYTHISANDPOS
 *		POLYMORPHICHROM
 *		POLYMORPHICEXP
 *	Invalid cases:
 */

TEST(uBasicNGSEXP_getOverlapping, VALIDEXP){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment* anotherExp = myExperiments.getExperiment("NoName_3elems");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*anotherExp); // TODO: Finish when the code is functional...
}

TEST(uBasicNGSEXP_getOverlapping, VALIDCHROM){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4"); 
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*aChrom); // TODO: Finish when the code is functional...
}

TEST(uBasicNGSEXP_getOverlapping, VALIDPOS){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping("chr4", 100, 200); // TODO: Finish when the code is functional...
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYEXP){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment* anotherExp = myExperiments.getExperiment("Empty_Exp");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*anotherExp); // TODO: Finish when the code is functional...
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYCHROM){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("NoName_Empty")->getpChrom(""); 
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*aChrom); // TODO: Finish when the code is functional...
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISEXP){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("Empty_Exp");
	uBasicNGSExperiment* anotherExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*anotherExp); // TODO: Finish when the code is functional...
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISCHROM){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("Empty_Exp");
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("MultipleChroms")->getpChrom(""); 
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*aChrom); // TODO: Finish when the code is functional...
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISANDPOS){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping("chr4", 1000, 2000); // TODO: Finish when the code is functional...
}

TEST(uBasicNGSEXP_getOverlapping, CHROMDONTEXISTS){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("Empty_Exp");
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("MultipleChroms")->getpChrom(""); 
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*aChrom); // TODO: Finish when the code is functional...
	EXPECT_EQ(overlapExp.count(), 0);
}

//TEST(uBasicNGSEXP_getOverlapping, POLYMORPHICHROM){ //TODO: What is a polymorphic chrom?
//}
//
//TEST(uBasicNGSEXP_getOverlapping, POLYMORPHICEXP){ //TODO: What is a polymorphic chrom?
//}

///**< set/get Chr Size */
//TEST(uBasicNGSEXP_getChrSize, VALID){
//       ASSERT_TRUE(false);
// }
// TEST(uBasicNGSEXP_getChrSize, NOTVALID){
//       ASSERT_TRUE(false);
// }
//  TEST(uBasicNGSEXP_getChrSize, EMPTY){
//       ASSERT_TRUE(false);
// }
//
// TEST(uBasicNGSEXP_setOverlapping, INVALIDVALUE){
//       ASSERT_TRUE(false);
// }
// TEST(uBasicNGSEXP_setOverlapping, INVALUECHR){
//       ASSERT_TRUE(false);
// }
// TEST(uBasicNGSEXP_setOverlapping, VALID){
//       ASSERT_TRUE(false);
// }
//
///**< getSubsetCount */
// TEST(uBasicNGSEXP_getSubsetCount, POSITIONS){
//       ASSERT_TRUE(false);
// }
//  TEST(uBasicNGSEXP_getSubsetCount, ELEMENT){
//       ASSERT_TRUE(false);
// }
// TEST(uBasicNGSEXP_getSubsetCount, CHROMNOEXIST){
//       ASSERT_TRUE(false);
// }
// TEST(uBasicNGSEXP_getSubsetCount, NONAMECHROM){
//       ASSERT_TRUE(false);
// }
