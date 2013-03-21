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
 * Tests for the function:
 *		InitialValue accumulateChromsInfo(BinaryOperation binary_op, InitialValue init) const;
 *	Valid cases:
 *		EMPTY
 *		NORMAL
 *	Invalid cases:
 *		EXCEPTION
 */
TEST(uBasicNGSEXP_STL_accumulateChromInfos, EMPTY){
	validExperiments myExperiments;
	int siteCount=0;
	auto functOp = [&](int siteCounts,const uBasicNGSChrom& item)
	{
		return (siteCounts+=item.sumSiteSize());
	};
	EXPECT_EQ(myExperiments.getExperiment("Empty_Exp")->accumulateChromsInfo(functOp,siteCount), 0);
}

TEST(uBasicEXP_STL_accumulateChromInfo, NORMAL){
	validExperiments myExperiments;
	int siteCount=0;
	auto functOp = [&](int siteCounts,const uBasicNGSChrom& item)
	{
		return (siteCounts+=item.sumSiteSize());
	};
	EXPECT_EQ(myExperiments.getExperiment("NoName_1elem")->accumulateChromsInfo(functOp,siteCount), 101);

}

TEST(uBasicEXP_STL_accumulateChromInfo, EXCEPTION){
	validExperiments myExperiments;
	auto functOp = [](int siteCounts,const uBasicNGSChrom & item)
	{
		throw param_throw();
		return 0;
	};
	int siteCount =0;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->accumulateChromsInfo(functOp,siteCount),param_throw);
}

/*
 * Tests for the function:
 *		auto computeOnAllChroms(UnaryOperation unary_op) const -> std::map<std::string, decltype(unary_op(_CHROM_()))>;
 *
 *	Valid cases:
 *		NORMAL
 *		EMPTY
 *	Invalid cases:
 *		THROW
 */

TEST(uBasicEXP_STL_computeOnAllChroms, NORMAL){
	validExperiments myExperiments;
	auto functOp = [&](uBasicNGSChrom item) -> unsigned long long
	{
		return item.sumSiteSize();
	};
	auto results = myExperiments.getExperiment("MultipleChroms")->computeOnAllChroms(functOp);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->sumSiteSize(), results.at(""));
}

TEST(uBasicEXP_STL_computeOnAllChroms, EMPTY){
	validExperiments myExperiments;
	auto functOp = [&](uBasicNGSChrom item) -> unsigned long long
	{
		return item.sumSiteSize();
	};
	auto results = myExperiments.getExperiment("Empty_Exp")->computeOnAllChroms(functOp);
	EXPECT_EQ((int)(results.size()), 0);
}

TEST(uBasicEXP_STL_computeOnAllChroms, THROW) {
	validExperiments myExperiments;
	auto functOp = [&](uBasicNGSChrom item) -> unsigned long long
	{
		throw param_throw();
		return 0;
	};
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->computeOnAllChroms(functOp),param_throw);
}

/*
 * Tests for the function:

 // TODO: How do you test applyOnAllChroms const??
 *		auto getSpecificChroms(UnaryPredicate pred) const->decltype(ExpMap)
 *	Valid Cases:
 *		NORMAL
 *		NOMATCH
 *		NOCHROM
 * 	Ivalid Cases:
 */
bool sumSiteSizeIs0(const uBasicNGSChrom & chrom)
{
	return chrom.sumSiteSize() == 0;
}

TEST(uBasicEXP_STL_getSpecificChroms, NORMAL) {
	validExperiments myExperiments;
	auto results = myExperiments.getExperiment("MultipleChroms")->getSpecificChroms(sumSiteSizeIs0);
	EXPECT_EQ((int)(results.size()), 1);
	EXPECT_EQ(results["chr3"].getChr(), myExperiments.getExperiment("MultipleChroms")->getpChrom("chr3")->getChr());
}

TEST(uBasicEXP_STL_getSpecificChroms, NOMATCH) {
	validExperiments myExperiments;
	auto results = myExperiments.getExperiment("NoName_1elem")->getSpecificChroms(sumSiteSizeIs0);
	EXPECT_EQ((int)(results.size()), 0);
}

TEST(uBasicEXP_STL_getSpecificChroms, NOCHROM) {
	validExperiments myExperiments;
	auto results = myExperiments.getExperiment("Empty_Exp")->getSpecificChroms(sumSiteSizeIs0);
	EXPECT_EQ((int)(results.size()), 0);
}

/*
 * Tests for the function:
 *		template<class UnaryFunction> UnaryFunction applyOnAllChroms(UnaryFunction f);
 *	Valid Cases:
 *		NORMAL
 *		NOCHROM
 * 	Ivalid Cases:
 *		THROW
 */

TEST(uBasicEXP_STL_applyOnAllChroms, NORMAL) {
	validExperiments myExperiments;
	auto functOp = [] (uBasicNGSChrom& chrom)
	{
		chrom.divideItemsIntoNBins(2, SplitType::IGNORE);
	};
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 13);
	myExperiments.getExperiment("MultipleChroms")->applyOnAllChroms(functOp);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 26);
}

TEST(uBasicEXP_STL_applyOnAllChroms, NOCHROM) {
	validExperiments myExperiments;
	auto functOp = [] (uBasicNGSChrom& chrom)
	{
		chrom.divideItemsIntoNBins(2, SplitType::IGNORE);
	};
	EXPECT_NO_THROW(myExperiments.getExperiment("Empty_Exp")->applyOnAllChroms(functOp));
}

TEST(uBasicEXP_STL_applyOnAllChroms, THROW) {
	validExperiments myExperiments;
	auto functOp = [] (uBasicNGSChrom& chrom)
	{
		throw param_throw();
	};
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->applyOnAllChroms(functOp), param_throw);
}

/*
 * Tests for the function:
 *		template<class UnaryFunction> UnaryFunction applyOnAllChroms(UnaryFunction f) const;
 *	Valid Cases:
 *		NORMAL
 *		NOCHROM
 * 	Ivalid Cases:
 *		THROW
 */

TEST(uBasicEXP_STL_applyOnAllChromsConst, NORMAL) {
	validExperiments myExperiments;
	int count = 0;
	auto functOp = [&](const uBasicNGSChrom& a)
	{
		count++;
	};
	EXPECT_EQ(count, 0);
	const uBasicNGSExperiment myConstExp = *(myExperiments.getExperiment("MultipleChroms"));
	myConstExp.applyOnAllChroms(functOp);
	EXPECT_EQ(count, 6);
}

TEST(uBasicEXP_STL_applyOnAllChromsConst, NOCHROM) {
	validExperiments myExperiments;
	int count = 0;
	auto functOp = [&](const uBasicNGSChrom& a)
	{
		count++;
	};
	EXPECT_EQ(count, 0);
	const uBasicNGSExperiment myConstExp = *(myExperiments.getExperiment("Empty_Exp"));
	myConstExp.applyOnAllChroms(functOp);
	EXPECT_EQ(count, 0);
}

TEST(uBasicEXP_STL_applyOnAllChromsConst, THROW) {
	validExperiments myExperiments;
	auto functOp = [&](const uBasicNGSChrom& a)
	{
		throw param_throw();
	};
	const uBasicNGSExperiment myConstExp = *(myExperiments.getExperiment("MultipleChroms"));
	EXPECT_THROW(myConstExp.applyOnAllChroms(functOp), param_throw);
}

/*
 * Tests for the function:
 *		template<class UnaryFunction>
 *		UnaryFunction applyOnOneChrom(UnaryFunction f, const std::string & chr);
 *
 *	Valid Cases:
 *		NORMAL
 *		NONAMECHROM
 * 	Ivalid Cases:
 *		THROW
 *		CHROMDONTEXISTS
 */
//
//TEST(uBasicEXP_STL_applyOnOneChrom, NORMAL) {
//	validExperiments myExperiments;
//	auto functOp = [](uBasicNGSChrom& item)
//	{
//		return item.divideItemsIntoNBins(2, SplitType::IGNORE);
//	};
//	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 13);
//	myExperiments.getExperiment("MultipleChroms")->applyOnOneChrom(functOp, "chr4");
//	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 14);
//}
//
//TEST(uBasicEXP_STL_applyOnOneChrom, NOCHROM) {
//	validExperiments myExperiments;
//	auto functOp = [](uBasicNGSChrom& item)
//	{
//		return item.divideItemsIntoNBins(2, SplitType::IGNORE);
//	};
//	EXPECT_THROW(myExperiments.getExperiment("Empty_Exp")->applyOnOneChrom(functOp, "chr4"), param_throw);
//
//}
//
//TEST(uBasicEXP_STL_applyOnOneChrom, THROW) {
//	validExperiments myExperiments;
//	auto functOp = [&](const uBasicNGSChrom& a)
//	{
//		throw param_throw();
//	};
//	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->applyOnOneChrom(functOp,"chr4"),param_throw);
//}
//
//TEST(uBasicEXP_STL_applyOnOneChrom, CHROMDONTEXISTS) {
//	validExperiments myExperiments;
//	int count = 0;
//	auto functOp = [&](const uBasicNGSChrom& a)
//	{
//		count++;
//	};
//	EXPECT_THROW(myExperiments.getExperiment("Empty_Exp")->applyOnOneChrom(functOp,"chr4"),param_throw);
//}

/*
 * Tests for the function:
 *		template<class UnaryFunction>
 *		UnaryFunction applyOnSites(UnaryFunction f);
 *
 *	Valid Cases:
 *		NORMAL
 *		EMPTY
 * 	Ivalid Cases:
 *		THROW
 */

TEST(uBasicEXP_STL_applyOnSites, NORMAL) {
	validExperiments myExperiments;
	auto functOp = [](uBasicNGS & item){ item.extendSite(20); };
	myExperiments.getExperiment("NoName_1elem")->applyOnSites(functOp);
	EXPECT_TRUE(myExperiments.getExperiment("NoName_1elem")->getpChrom("")->getSite(0).isEqual(uBasicNGS("",80,220)));
}

TEST(uBasicEXP_STL_applyOnSites, EMPTY) {
	validExperiments myExperiments;
	auto functOp = [](uBasicNGS & item){ item.extendSite(20); };
	EXPECT_NO_THROW(myExperiments.getExperiment("Empty_Exp")->applyOnSites(functOp));
}

TEST(uBasicEXP_STL_applyOnSites, THROW) {
	validExperiments myExperiments;
	auto functOp = [](uBasicNGS & item) { throw param_throw(); };
	EXPECT_THROW(myExperiments.getExperiment("NoName_1elem")->applyOnSites(functOp), param_throw);
}

/*
 * Tests for the function:
 *		template<class UnaryFunction>
 *		UnaryFunction applyOnSites(const UnaryFunction f)const;
 *
 *	Valid Cases:
 *		NORMAL
 *		EMPTY
 * 	Ivalid Cases:
 *		THROW
 */

TEST(uBasicEXP_STL_applyOnSitesConst, NORMAL) {
	validExperiments myExperiments;
	int count = 0;
	auto functOp = [&](const uBasicNGS& a)
	{
		count++;
	};
//	auto functOp = [&](uBasicNGS item)->int{  return item.getLenght();};
	EXPECT_EQ(count, 0);
//	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->applyOnSites(functOp));
	const uBasicNGSExperiment myConstExp = *(myExperiments.getExperiment("MultipleChroms"));
	myConstExp.applyOnSites(functOp); // TODO: Weird warning
	EXPECT_EQ(count, 13);
}

/*
 * Tests for the function:
 *		template<class UnaryFunction>
 *		void loadWithParserAndRun(std::string filepath, std::string pType, UnaryFunction f, int pBlockSize=1);
 *
 *	Valid Cases:
 *		VALIDSTREAM
 *	Invalid Cases:
 *
 */

TEST(uBasicEXP_STL_loadWithParserAndRunParser, VALIDSTREAM) {
	uParser aParser("../data/BED/test.bed", "BED");
	auto functOp = [&](uBasicNGS  item)
	{
	};
	uBasicNGSExperiment anExp;
	anExp.loadWithParserAndRun(aParser, functOp);
}

/*
 * Tests for the function:
 *		countChromsWithProperty(UnaryPredicate pred) const;
 *		template<class Compare>
 *
 *	Valid Cases:
 *	Invalid Cases:
 *		THROW
 */
/*
 * Tests for the function:
 *		NGSExpConstIter maxChrom(Compare comp) const;
 *		template<class Compare>
 *
 *	Valid Cases:
 *	Invalid Cases:
 *		THROW
 */
/*
 * Tests for the function:
 *		NGSExpConstIter minChrom(Compare comp) const;
 *		template<class Compare>
 *
 *	Valid Cases:
 *	Invalid Cases:
 *		THROW
 */
/*
 * Tests for the function:
 *		std::pair<NGSExpConstIter, NGSExpConstIter> minAndMaxChroms(Compare comp) const;
 *
 *	Valid Cases:
 *	Invalid Cases:
 *		THROW
 */
