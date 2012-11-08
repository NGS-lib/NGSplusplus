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

