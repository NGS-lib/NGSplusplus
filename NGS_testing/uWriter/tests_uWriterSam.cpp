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
TEST_F(TestSamWriter, WriteToken_ValidSam) {
	string expected = "ab00001\t256	chr1\t1\t255\t2M3I2X1=12X\t*\t*\t0\tACGTN.acgtn.ACGTGTCN\t####################";

   // writerSAM.writeHeader();


	ASSERT_TRUE(m_pOssSAM->str().find(expected) != string::npos);
}

TEST(TestSamWriter_NoF, WriteValidHeader) {
    /**< Initialize streams */
   //  std::ostream & outFile = ((outputPath.size()!=0) ? outputOS : std::cout);

    auto m_pOssSAM = new ostringstream(ostringstream::out);
    /**< Initialize writers */
    uWriter writerSAM(m_pOssSAM, "SAM");

    writerSAM.addToHeader(header_param::CHR,"chr1");
    writerSAM.addToHeader(header_param::CHR_SIZE,"200513");

    writerSAM.writeHeader();

    cout <<m_pOssSAM->str();
//	ASSERT_TRUE(m_pOssSAM->str().find(expected) != string::npos);
}

