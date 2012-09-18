#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "IO/Writer/uWriter.h"
#include "tests_Writer_Fixtures.h"

using namespace std;
using namespace NGS;

/*
 * Tests for the function: 
 *		void writeToken(const uToken& token);
 *	Valid case:
 *		SixValidValues
 *		NoScoreValue
 *		NoNameValue
 *		NoScoreStrandValues
 */
TEST_F(TestsBedWriter, WriteToken_SixValidValues) {
	string expected = "chr1\t1\t21\tab00001\t111\t+\n";
	ASSERT_TRUE(m_pOss->str().find(expected) != string::npos);
}

TEST_F(TestsBedWriter, WriteToken_NoScoreValue) {
	string expected = "chr1\t1\t21\tab00001\t.\t+\n";
	ASSERT_TRUE(m_pOss->str().find(expected) != string::npos);
}

TEST_F(TestsBedWriter, WriteToken_NoNameValue) {
	string expected = "chr1\t1\t21\t.\t111\t+\n";
	ASSERT_TRUE(m_pOss->str().find(expected) != string::npos);
}

TEST_F(TestsBedWriter, WriteToken_NoScoreStrandValues) {
	string expected = "chr1\t1\t21\tab00001\n";
	ASSERT_TRUE(m_pOss->str().find(expected) != string::npos);
}

