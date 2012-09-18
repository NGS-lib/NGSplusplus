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
 *		NormalUsage
 */
TEST_F(TestsBedWriter, WriteToken_NormalUsage) {
	string expected = "chr1\t1\t21\tab00001\t.\t+\n";
	expected.append("chr2\t101\t121\tab00002\t.\t-\n");
	ASSERT_EQ(m_pOss->str(), expected);
}
