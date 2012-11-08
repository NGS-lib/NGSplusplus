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
 *		AllValuesArePresent
 *		MissingValidValue
 *		OneInvalidValue
 *		NoValidValues
 */

TEST_F(TestsCustomWriter, WriteToken_AllValuesArePresent) {
	string expected = "chr1\t1\tACGTN.acgtn.ACGTGTCN\n";
	ASSERT_TRUE(m_pOss->str().find(expected) != string::npos);
}

TEST_F(TestsCustomWriter, WriteToken_MissingValidValue) {
	string expected = "chr1\t1\t.\n";
	ASSERT_TRUE(m_pOss->str().find(expected) != string::npos);
}

TEST_F(TestsCustomWriterOneInvalidValue, WriteToken_OneInvalidValue) {
	string expected = "chr1\t1\tACGTN.acgtn.ACGTGTCN\t.\n";
	ASSERT_TRUE(m_pOss->str().find(expected) != string::npos);
}

TEST_F(TestsCustomWriterNoValidValues, WriteToken_NoValidValues) {
	string expected = ".\t.\t.\n";
	ASSERT_TRUE(m_pOss->str().find(expected) != string::npos);
}
