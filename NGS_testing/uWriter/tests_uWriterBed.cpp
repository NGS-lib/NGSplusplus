#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "NGS++.h"
#include "tests_Writer_Fixtures.h"

using namespace std;
using namespace NGS;
/*
 * Tests for the function:
 *		void writeToken(const uToken& token);
 *	Valid case:
 *		FourValidValuesBed4
 *		SixValidValuesBed6
 *		NoScoreValueBed6
 *		NoNameValueBed4
 *		NoNameValueBed6
 *		NoScoreStrandValuesBed6
 */
TEST_F(TestsBedWriter, WriteToken_FourValidValuesBed4) {
	string expected = "chr1\t1\t21\tab00001\n";
	//cout << m_pOssBed4->str();
	ASSERT_TRUE(m_pOssBed4->str().find(expected) != string::npos);
}

TEST_F(TestsBedWriter, WriteToken_SixValidValuesBed6) {
	string expected = "chr1\t1\t21\tab00001\t111\t+\n";
	ASSERT_TRUE(m_pOssBed6->str().find(expected) != string::npos);
}
TEST_F(TestsBedWriter, WriteToken_NoScoreValueBed6) {
	string expected = "chr1\t1\t21\tab00001\t.\t+\n";
	ASSERT_TRUE(m_pOssBed6->str().find(expected) != string::npos);
}

TEST_F(TestsBedWriter, WriteToken_NoNameValueBed4) {
	string expected = "chr1\t1\t21\t.\n";
	ASSERT_TRUE(m_pOssBed4->str().find(expected) != string::npos);
}

TEST_F(TestsBedWriter, WriteToken_NoNameValueBed6) {
	string expected = "chr1\t1\t21\t.\t111\t+\n";
	ASSERT_TRUE(m_pOssBed6->str().find(expected) != string::npos);
}
TEST_F(TestsBedWriter, WriteToken_NoScoreStrandValuesBed6) {
	string expected = "chr1\t1\t21\tab00001\t.\t.\n";
	ASSERT_TRUE(m_pOssBed6->str().find(expected) != string::npos);
}



ostringstream* m_pOsBed = new ostringstream(ostringstream::out);
std::string bedExpect ="chr1\t100\t200\t.\n";

TEST(TestsBedWriter_WriteToken,TAG) {

	uTags testTag("chr1", 100, 200);
	uWriter bedWriter(m_pOsBed,"BED");
    testTag.writeToOutput(bedWriter);
	ASSERT_TRUE(m_pOsBed->str().find(bedExpect) != string::npos);
}

TEST(TestsBedWriter_WriteToken,REGION) {

	uRegion testRegion("chr1", 100, 200);
	uWriter bedWriter(m_pOsBed,"BED");
    testRegion.writeToOutput(bedWriter);
	ASSERT_TRUE(m_pOsBed->str().find(bedExpect) != string::npos);
}

TEST(TestsBedWriter_WriteToken,BASICNGS) {

	uBasicNGS testBasic("chr1", 100, 200);
	uWriter bedWriter(m_pOsBed,"BED");
	 testBasic.writeToOutput(bedWriter);
	ASSERT_TRUE(m_pOsBed->str().find(bedExpect) != string::npos);
}

TEST(TestsBedWriter_WriteToken,GENE) {

	uGene testUgene("chr1", 100, 200);
	uWriter bedWriter(m_pOsBed,"BED");
	testUgene.writeToOutput(bedWriter);

	ASSERT_TRUE(m_pOsBed->str().find(bedExpect) != string::npos);
}

