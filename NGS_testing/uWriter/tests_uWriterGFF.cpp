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
 */
TEST_F(TestGFFWriter, WriteGFFCol2) {
	string expected = "ab00001	.	.	1	21	.	+	.";
    ASSERT_TRUE(m_pOssGFF->str().find(expected) != string::npos);
}
TEST_F(TestGFFWriter, WriteGFFCol3) {
	string expected = "chr1	.	.	1	21	111	+	.";
    ASSERT_TRUE(m_pOssGFF->str().find(expected) != string::npos);
}
