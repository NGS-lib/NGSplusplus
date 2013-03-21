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
const string HG18UCSC="../data/GENEPRED/hg18GenePred.genePred";

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
	EXPECT_THROW(myExperiments.getChrom(""), param_throw);
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
	EXPECT_THROW(anExp.getpChrom(""), param_throw);
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
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->getSite("chr4", 1), param_throw);
}

TEST(uBasicNGSEXP_getSite, BELOW0){
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->getSite("chr4", -1), param_throw);
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
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*anotherExp);
	EXPECT_EQ(overlapExp.count(), 3);
}

TEST(uBasicNGSEXP_getOverlapping, VALIDCHROM){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*aChrom);
	EXPECT_EQ(overlapExp.count(), 1);
}

TEST(uBasicNGSEXP_getOverlapping, VALIDPOS){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping("chr4", 100, 200);
	EXPECT_EQ(overlapExp.count(), 1);
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYEXP){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment* anotherExp = myExperiments.getExperiment("Empty_Exp");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*anotherExp);
	EXPECT_EQ(overlapExp.count(), 0);
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYCHROM){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("NoName_Empty")->getpChrom("");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*aChrom);
	EXPECT_EQ(overlapExp.count(), 0);
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISEXP){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("Empty_Exp");
	uBasicNGSExperiment* anotherExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*anotherExp);
	EXPECT_EQ(overlapExp.count(), 0);
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISCHROM){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("Empty_Exp");
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("MultipleChroms")->getpChrom("");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*aChrom);
	EXPECT_EQ(overlapExp.count(), 0);
}

TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISANDPOS){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("MultipleChroms");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping("chr4", 1000, 2000);
	EXPECT_EQ(overlapExp.count(), 0);
}

TEST(uBasicNGSEXP_getOverlapping, CHROMDONTEXISTS){
	validExperiments myExperiments;
	uBasicNGSExperiment* anExp = myExperiments.getExperiment("Empty_Exp");
	uBasicNGSChrom* aChrom = myExperiments.getExperiment("MultipleChroms")->getpChrom("");
	uBasicNGSExperiment overlapExp = anExp->getOverlapping(*aChrom);
	EXPECT_EQ(overlapExp.count(), 0);
}

//TEST(uBasicNGSEXP_getOverlapping, POLYMORPHICEXP){
//	validExperiments myExperiments;
//	uGeneExperiment anExp;
//	uParser genePredParser(HG18UCSC,"GENEPRED");
//	anExp.loadWithParser(genePredParser,1);
//	uBasicNGSExperiment* anotherExp = myExperiments.getExperiment("NoName_3elems");
//	uBasicNGSExperiment overlapExp = anExp.getOverlapping(*anotherExp);
//	EXPECT_EQ(overlapExp.count(), 1);
//}


/*
 * Test for the function:
 *		int getChrSize(std::string chr);
 *	Valid cases:
 *		VALIDCHROM
 *		EMPTYCHROM
 *	Invalid cases:
 *		INVALIDCHROM
 */

TEST(uBasicNGSEXP_getChrSize, VALID){
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize("chr4"), 0);
	myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->inferChrSize();
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize("chr4"), 200);
}

TEST(uBasicNGSEXP_getChrSize, EMPTYCHROM){
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("NoName_Empty")->getChrSize(""), 0);
}

TEST(uBasicNGSEXP_getChrSize, INVALIDCHROM){
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("NoName_Empty")->getChrSize("chr4"), param_throw);
}

/*
 * Test for the functions:
 *		void setChrSize(std::string chr, int chrSize);
 *	Valid cases:
 *		VALIDCHRSIZE
 *		CHRSIZEZERO
 *		NONAMECHR
 *	Invalid cases:
 *		INVALIDCHRSIZE
 *		INVALIDCHRNAME
 */

TEST(uBasicNGSEXP_setChrSize, VALIDCHRSIZE){
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize("chr4"), 0);
	myExperiments.getExperiment("MultipleChroms")->setChrSize("chr4", 300);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize("chr4"), 300);
}

TEST(uBasicNGSEXP_setChrSize, CHRSIZEZERO){
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->setChrSize("chr4", 0);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize("chr4"), 0);
}

TEST(uBasicNGSEXP_setChrSize, NONAMECHR){
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->setChrSize("", 300);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize(""), 300);
}

TEST(uBasicNGSEXP_setChrSize, INVALIDCHRSIZE){
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->setChrSize("chr4", -1), param_throw);
}

TEST(uBasicNGSEXP_setChrSize, INVALIDCHRNAME){
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("NoName_Empty")->setChrSize("chr4", 200), param_throw);
}

/*
 * Tests for the function:
 *	     _CHROM_ getSubset(std::string pChr, float pStart, float pEnd, OverlapType options=OverlapType::OVERLAP_PARTIAL);
 *
 *	Valid cases:
 *		VALIDCHROMNAME
 *		NONAMECHROM
 *		NOELEMINSUBSET
 *		CUSTOMSORT
 *	Invalid cases:
 *		CHROMDONTEXISTS
 *		NOTSORTED
 */

TEST(uBasicNGSEXP_getSubset, VALIDCHROMNAME) {
	validExperiments myExperiments;
//	myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->sortSites();
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	uBasicNGSChrom aChrom = myExperiments.getExperiment("MultipleChroms")->getSubset("chr4", 100, 200);
	EXPECT_EQ(aChrom.count(), 1);
	aChrom.sortSites();
	auto it = aChrom.findNextSite(0);
	EXPECT_EQ(it->getChr(), "chr4");
	EXPECT_EQ(it->getStart(), 100);
	EXPECT_EQ(it->getEnd(), 200);
	EXPECT_EQ(aChrom.getSite(0).getChr(), "chr4");
	EXPECT_EQ(aChrom.getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom.getSite(0).getEnd(), 200);
}

TEST(uBasicNGSEXP_getSubset, NONAMECHROM) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->getpChrom("")->sortSites();
	uBasicNGSChrom aChrom = myExperiments.getExperiment("MultipleChroms")->getSubset("", 150, 250);

	EXPECT_EQ(aChrom.getSite(0).getChr(), "");
	EXPECT_EQ(aChrom.getSite(0).getStart(), 100);
	EXPECT_EQ(aChrom.getSite(0).getEnd(), 200);
	EXPECT_EQ(aChrom.getSite(1).getChr(), "");
	EXPECT_EQ(aChrom.getSite(1).getStart(), 200);
	EXPECT_EQ(aChrom.getSite(1).getEnd(), 300);
	EXPECT_THROW(aChrom.getSite(2), param_throw);
}

TEST(uBasicNGSEXP_getSubset, NOELEMINSUBSET) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->getpChrom("")->sortSites();
	uBasicNGSChrom aChrom = myExperiments.getExperiment("MultipleChroms")->getSubset("", 1000, 2000);

	EXPECT_THROW(aChrom.getSite(0), param_throw);
}

TEST(uBasicNGSEXP_getSubset, CUSTOMSORT) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr6")->compareLenght, &uBasicNGS::getLenght);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr6")->count(), 3);
	uBasicNGSChrom aChrom = myExperiments.getExperiment("MultipleChroms")->getSubset("chr6", 120, 170);
	EXPECT_EQ(aChrom.count(), 1);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr6")->count(), 3);
}

TEST(uBasicNGSEXP_getSubset, CHROMDONTEXISTS) {
	validExperiments myExperiments;
	uBasicNGSChrom aChrom = myExperiments.getExperiment("NoName_Empty")->getSubset("chr1", 1000, 2000);
	EXPECT_THROW(aChrom.getSite(0), param_throw);
}

TEST(uBasicNGSEXP_getSubset, NOTSORTED) {
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->getSubset("", 1000, 2000), unsorted_throw);
}

/*
 * Tests for the function:
 *		_SELF_ getDistinct( std::string pChr, float pStart, float pEnd, OverlapType type=OverlapType::OVERLAP_PARTIAL);
 *	Valid Cases:
 *		VALIDCHROMNAME
 *		NONAMECHROM
 *		NOELEMINSUBSET
 *		CUSTOMSORT
 *	Invalid Cases:
 *		CHROMDONTEXISTS
 *		NOTSORTED
 */

TEST(uBasicNGSEXP_getDistinct, VALIDCHROMNAME) {
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 13);
	myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->sortSites();
	uBasicNGSExperiment distinctExp = myExperiments.getExperiment("MultipleChroms")->getDistinct("chr4", 100, 200);
	EXPECT_EQ(distinctExp.count(), 12);
	EXPECT_THROW(distinctExp.getpChrom("chr4")->getSite(0), param_throw);
}

TEST(uBasicNGSEXP_getDistinct, NONAMECHROM) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->getpChrom("")->sortSites();
	uBasicNGSExperiment distinctExp = myExperiments.getExperiment("MultipleChroms")->getDistinct("", 100, 250);
	EXPECT_EQ(distinctExp.count(), 11);
	EXPECT_NO_THROW(distinctExp.getpChrom(""));
}

TEST(uBasicNGSEXP_getDistinct, CUSTOMSORT) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr6")->compareLenght, &uBasicNGS::getLenght);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr6")->count(), 3);
	uBasicNGSExperiment distinctExp = myExperiments.getExperiment("MultipleChroms")->getDistinct("chr6", 120, 170);
	EXPECT_EQ(distinctExp.getpChrom("chr6")->count(), 2);
	EXPECT_EQ(distinctExp.getpChrom("")->count(), 3);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr6")->count(), 3);
}

TEST(uBasicNGSEXP_getDistinct, NOELEMINSUBSET) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->getpChrom("")->sortSites();
	uBasicNGSExperiment distinctExp = myExperiments.getExperiment("MultipleChroms")->getDistinct("", 1000, 2000);
	EXPECT_EQ(distinctExp.count(), 13);
}

TEST(uBasicNGSEXP_getDistinct, CHROMDONTEXISTS) {
	validExperiments myExperiments;
	uBasicNGSExperiment distinctExp = myExperiments.getExperiment("NoName_Empty")->getDistinct("", 1000, 2000);
	EXPECT_THROW(distinctExp.getpChrom("chr1"), param_throw);
}

TEST(uBasicNGSEXP_getDistinct, NOTSORTED) {
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->getDistinct("", 1000, 2000), unsorted_throw);
}

/*
 * Test for the function:
 *		int getSubsetCount(const std::string & chr, const float start, const float end, const OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
 *		int getSubsetCount(const _BASE_ & subsetReg, const OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
 *	Valid cases:
 *		POSITIONS
 *		ELEMENT
 *		NONAMECHROM
 *	Invalid cases:
 *		CHROMNOEXIST
 *		NOTSORTED //TODO: unsorted_throw
 */

TEST(uBasicNGSEXP_getSubsetCount, POSITIONS){
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->sortSites();
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSubsetCount("chr4", 100, 200), 1);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSubsetCount("chr4", 1000, 2000), 0);
}

TEST(uBasicNGSEXP_getSubsetCount, ELEMENT){
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->sortSites();
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSubsetCount(uBasicNGS("chr4",100,200)), 1);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSubsetCount(uBasicNGS("chr4",1000,2000)), 0);
}

TEST(uBasicNGSEXP_getSubsetCount, NONAMECHROM){
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->getpChrom("")->sortSites();
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getSubsetCount("", 100, 200), 2);
}

TEST(uBasicNGSEXP_getSubsetCount, CHROMNOEXIST){
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->sortSites();
	EXPECT_THROW(myExperiments.getExperiment("NoName_Empty")->getSubsetCount("chr4", 100, 200), param_throw);
}

/*
 * Tests for the function:
 *		void setChrSize(std::string chr, int chrSize);
 *		int getChrSize(std::string chr);
 *	Valid cases:
 *		CHRSIZEWASNOTSET
 *		CHRSIZEWASSET
 *	Invalid cases:
 *		CHROMDONTEXISTS
 */

TEST(uBasicNGSEXP_setChrSize, CHRSIZEWASNOTSET) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->setChrSize("", 2000);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize(""), 2000);
}

TEST(uBasicNGSEXP_setChrSize, CHRSIZEWASSET) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->setChrSize("", 2000);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize(""), 2000);
	myExperiments.getExperiment("MultipleChroms")->setChrSize("", 200);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getChrSize(""), 200);
}

TEST(uBasicNGSEXP_setChrSize, CHROMDONTEXISTS) {
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("NoName_Empty")->setChrSize("chr1", 2000), param_throw);
}

/*
 * Tests for the function:
 *		void divideItemsIntoBinofSize(int N, SplitType type=SplitType::STRICT);
 *	Valid cases:
 *		EXPWITHELEMS
 *		EXPWITHOUTELEM
 *	Invalid cases:
 */

TEST(uBasicNGSEXP_divideItemsIntoBinofSize, EXPWITHELEMS) {
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("NoName_1elem")->count(), 1);
	myExperiments.getExperiment("NoName_1elem")->divideItemsIntoBinofSize(51, SplitType::ADD);
	EXPECT_EQ(myExperiments.getExperiment("NoName_1elem")->count(), 2);
}

TEST(uBasicNGSEXP_divideItemsIntoBinofSize, EXPWITHOUTELEM) {
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("NoName_Empty")->count(), 0);
	myExperiments.getExperiment("NoName_Empty")->divideItemsIntoBinofSize(51, SplitType::ADD);
	EXPECT_EQ(myExperiments.getExperiment("NoName_Empty")->count(), 0);
}

/*
 * Tests for the function:
 *		void divideItemsIntoNBins(int N, SplitType type=SplitType::STRICT);
 *	Valid cases:
 *		EXPWITHELEMS
 *		EXPWITHOUTELEM
 *	Invalid cases:
 */

TEST(uBasicNGSEXP_divideItemsIntoNBins, EXPWITHELEMS) {
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 13);
	myExperiments.getExperiment("MultipleChroms")->divideItemsIntoNBins(2, SplitType::IGNORE);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 26);
}

TEST(uBasicNGSEXP_divideItemsIntoNBins, EXPWITHOUTELEM) {
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("NoName_Empty")->count(), 0);
	myExperiments.getExperiment("NoName_Empty")->divideItemsIntoNBins(51, SplitType::ADD);
	EXPECT_EQ(myExperiments.getExperiment("NoName_Empty")->count(), 0);
}

/*
 * Tests for the function:
 *		template<class UnaryPredicate>
 *		void removeSpecificSites(UnaryPredicate pred);
 *
 * 	Valid cases:
 *		NONREMOVED
 *		SOMEREMOVED
 *		ALLREMOVED
 *		EMPTY
 *	Invalid cases:
 */
TEST(uBasicNGSGENEXP_removeSpecificSites, NONREMOVED) {
	validExperiments myExperiments;
	auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>2000); };
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 13);
	myExperiments.getExperiment("MultipleChroms")->removeSpecificSites(functOp);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 13);
}

TEST(uBasicNGSGENEXP_removeSpecificSites, SOMEREMOVED) {
	validExperiments myExperiments;
	auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>101); };
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 13);
	myExperiments.getExperiment("MultipleChroms")->removeSpecificSites(functOp);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 10);
}

TEST(uBasicNGSGENEXP_removeSpecificSites, ALLREMOVED) {
	validExperiments myExperiments;
	auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>99); };
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 13);
	myExperiments.getExperiment("MultipleChroms")->removeSpecificSites(functOp);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->count(), 0);
}

TEST(uBasicNGSGENEXP_removeSpecificSites, EMPTY) {
	validExperiments myExperiments;
	auto functOp = [](const uBasicNGS & item){  return (item.getLenght()>99); };
	EXPECT_EQ(myExperiments.getExperiment("Empty_Exp")->count(), 0);
	myExperiments.getExperiment("Empty_Exp")->removeSpecificSites(functOp);
	EXPECT_EQ(myExperiments.getExperiment("Empty_Exp")->count(), 0);
}

/*
 * Tests for the function:
 *		void replaceChr(const _CHROM_ &);
 *
 * 	Valid cases:
 *		CHREXISTS
 *		CHRDOESNTEXISTS
 *	Invalid cases:
 */

TEST(uBasicNGSGENEXP_replaceChr, CHREXISTS) {
	uBasicNGSChrom aChr;
	aChr.setChr("");
	aChr.addData(uBasicNGS("",100,200));
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->count(), 3);
	myExperiments.getExperiment("MultipleChroms")->replaceChr(aChr);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->count(), 1);
}

TEST(uBasicNGSGENEXP_replaceChr, CHRDOESNTEXISTS) {
	uBasicNGSChrom aChr;
	aChr.setChr("chrY");
	aChr.addData(uBasicNGS("chrY",100,200));
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->replaceChr(aChr);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chrY")->count(), 1);
}

/*
 * Tests for the function:
 *		void removeChr(const std::string &);
 *
 * 	Valid cases:
 *		CHREXISTS
 *		CHRDOESNTEXISTS
 *	Invalid cases:
 */

TEST(uBasicNGSGENEXP_removeChr, CHREXISTS) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->removeChr("");
	EXPECT_FALSE(myExperiments.getExperiment("MultipleChroms")->isChrom(""));
}

TEST(uBasicNGSGENEXP_removeChr, CHRDOESNTEXISTS) {
	validExperiments myExperiments;
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->removeChr("chrY"));
}

/*
 * Tests for the function:
 *		_CHROM_ removeSubset(const std::string & pChr,const double pStart, const double pEnd, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
 *
 * 	Valid cases:
 *		CHREXISTS
 *		CHRDOESNTEXISTS
 *	Invalid cases:
 *		NOTSORTED
 */

TEST(uBasicNGSGENEXP_removeSubset, CHREXISTS) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->count(), 1);
	myExperiments.getExperiment("MultipleChroms")->removeSubset("chr4", 150, 151);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->count(), 0);
}

TEST(uBasicNGSGENEXP_removeSubset, CHRDOESNTEXISTS) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->removeSubset("chrY", 150, 151));
}

TEST(uBasicNGSGENEXP_removeSubset, NOTSORTED) {
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->removeSubset("chr4", 150, 151), unsorted_throw);
}

/*
 * Tests for the function:
 *		_SELF_ removeDistinct(const std::string & pChr,const double pStart, const double pEnd, OverlapType options=OverlapType::OVERLAP_PARTIAL);
 *
 * 	Valid cases:
 *		CHREXISTS
 *		CHRDOESNTEXISTS
 *	Invalid cases:
 *		NOTSORTED
 */

TEST(uBasicNGSGENEXP_removeDistinct, CHREXISTS) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->count(), 1);
	auto distinctExp = myExperiments.getExperiment("MultipleChroms")->removeDistinct("chr4", 150, 151);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->count(), 1);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->isChrom("chr4"), true);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->isChrom(""), false);
	EXPECT_EQ(distinctExp.isChrom("chr4"), true);
	EXPECT_EQ(distinctExp.getpChrom("chr4")->count(), 0);
}

TEST(uBasicNGSGENEXP_removeDistinct, CHROMDONTEXISTS) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	auto distinctExp = myExperiments.getExperiment("MultipleChroms")->removeDistinct("chrY", 150, 151);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->isChrom("chr4"), false);
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->isChrom(""), false);
	EXPECT_EQ(distinctExp.isChrom("chrY"), false);
}

TEST(uBasicNGSGENEXP_removeDistinct, NOTSORTED) {
	validExperiments myExperiments;
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->removeDistinct("chr4", 150, 151), unsorted_throw);
}

/*
 * Tests for the function:
 *		virtual void addData(const _BASE_ &);
 * 	Valid cases:
 *		VALID
 *	Invalid cases:
 */
TEST(uBasicNGSGENEXP_addDataBase, VALID) {
	validExperiments myExperiments;
	uBasicNGS base("chr1", 100, 200);
	EXPECT_NO_THROW(myExperiments.getExperiment("Empty_Exp")->addData(base));
	EXPECT_EQ(myExperiments.getExperiment("Empty_Exp")->getpChrom("chr1")->count(), 1);
}

/*
 * Tests for the function:
 *		virtual void addData(const _CHROM_ &);
 * 	Valid cases:
 *		CHREXISTS
 *		CHROMDONTEXISTS
 *	Invalid cases:
 */

TEST(uBasicNGSGENEXP_addDataChrom, CHREXISTS) {
	validExperiments myExperiments;
	uBasicNGSChrom aChr;
	aChr.setChr("chr3");
	aChr.addData(uBasicNGS("chr3",100,200));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr3")->count(), 0);
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->addData(aChr));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr3")->count(), 1);
}

TEST(uBasicNGSGENEXP_addDataChrom, CHRDOESNTEXISTS) {
	validExperiments myExperiments;
	uBasicNGSChrom aChr;
	aChr.setChr("chr1");
	aChr.addData(uBasicNGS("chr1",100,200));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->isChrom("chr1"), false);
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->addData(aChr));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr1")->count(), 1);
}
/*
 * Tests for the function:
 *		virtual void addData(const uToken &);
 * 	Valid cases:
 *		VALID
 *	Invalid cases:
 */

TEST(uBasicNGSGENEXP_addDataToken, VALID) {
	uBasicNGS uTest("chr3", 1000, 1200);
	uToken ourToken= uTest.createToken();
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr3")->count(), 0);
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->addData(ourToken));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr3")->count(), 1);
}

/*
 * Tests for the function:
 *		void removeSite(const std::string & pChr,const long int position);
 * 	Valid cases:
 *		VALID
 *	Invalid cases:
 *		INDEXTOOHIGH
 *		EMPTYCHR
 */
/* // Removed since they are now protected
TEST(uBasicNGSGENEXP_removeSitePosition, VALID) {
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->count(), 1);
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->removeSite("chr4", 0));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->count(), 0);
}

TEST(uBasicNGSGENEXP_removeSitePosition, INDEXTOOHIGH) {
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr4")->count(), 1);
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->removeSite("chr4", 1), param_throw);
}

TEST(uBasicNGSGENEXP_removeSitePosition, INDEXTOOHIGH) {
	validExperiments myExperiments;
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("chr3")->count(), 0);
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->removeSite("chr3", 0), param_throw);
}
*/
/*
 * Tests for the function:
 *		void removeSite(const std::string & pChr,const long int pStart,const long int pEnd);
 * 	Valid cases:
 *		VALID
 *	Invalid cases:
 */

// Protected

//TEST(uBasicNGSGENEXP_removeSiteBeginEnd, VALID) {
//	ASSERT_TRUE(false);
//}

/*
 * Tests for the function:
 *		void removeSite(const std::string & pChr,VecGenConstIter pItrPos);
 * 	Valid cases:
 *		VALID
 *	Invalid cases:
 */

TEST(uBasicNGSGENEXP_removeSiteIter, VALID) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	auto itr = myExperiments.getExperiment("MultipleChroms")->findNextSite("", 1);
	EXPECT_EQ(itr->getChr(), "");
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->count(), 3);
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->removeSite(itr));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->count(), 2);
}

/*
 * Tests for the function:
 *		void removeSite(const std::string & pChr,VecGenConstIter pItrStart,VecGenConstIter pItrEnd);
 * 	Valid cases:
 *		VALID
 *		SAMESITE
 *	Invalid cases:
 *		ITERNOTONSAMECHROM
 *		NOTSORTED
 */

TEST(ubasicNGSGENEEXP_removeSiteIterBeginEnd, VALID) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	auto itr1 = myExperiments.getExperiment("MultipleChroms")->findNextSite("", 1);
	auto itr2 = myExperiments.getExperiment("MultipleChroms")->findNextSite("", 150);
	EXPECT_EQ(itr1->getChr(), "");
	EXPECT_EQ(itr2->getChr(), "");
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->count(), 3);
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->removeSite(itr1, itr2));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->count(), 2);
}

TEST(ubasicNGSGENEEXP_removeSiteIterBeginEnd, SAMESITE) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	auto itr1 = myExperiments.getExperiment("MultipleChroms")->findNextSite("", 1);
	auto itr2 = myExperiments.getExperiment("MultipleChroms")->findNextSite("", 1);
	EXPECT_EQ(itr1->getChr(), "");
	EXPECT_EQ(itr2->getChr(), "");
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->count(), 3);
	EXPECT_NO_THROW(myExperiments.getExperiment("MultipleChroms")->removeSite(itr1, itr2));
	EXPECT_EQ(myExperiments.getExperiment("MultipleChroms")->getpChrom("")->count(), 3);
}

TEST(uBasicNGSGENEXP_removeSiteIterBeginEnd, ITERNOTONSAMECHROM) {
	validExperiments myExperiments;
	myExperiments.getExperiment("MultipleChroms")->sortSites();
	auto itr1 = myExperiments.getExperiment("MultipleChroms")->findNextSite("", 1);
	auto itr2 = myExperiments.getExperiment("MultipleChroms")->findNextSite("chr4", 1);
	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->removeSite(itr1, itr2), param_throw);
}

// NOTE: Cannot test unsorted_throw on removeSite because it throws when I try to get an iterator with findNextSite on unsorted exp.
//TEST(uBasicNGSGENEXP_removeSiteIterBeginEnd, NOTSORTED) {
//	validExperiments myExperiments;
//	auto itr1 = myExperiments.getExperiment("MultipleChroms")->findNextSite("", 1);
//	auto itr2 = myExperiments.getExperiment("MultipleChroms")->findNextSite("", 1);
//	EXPECT_THROW(myExperiments.getExperiment("MultipleChroms")->removeSite(itr1, itr2), unsorted_throw);
//}

// Tested indirectly
//TEST(uBasicNGSGENEXP_sortSites, VALID) {
//	ASSERT_TRUE(false);
//}
