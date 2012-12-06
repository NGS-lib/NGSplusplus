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

TEST(TestGFFWriterConvert, FromBedReadGFF) {
   //  auto p_SSGFF = new ostringstream(ostringstream::out);
    uParser Parser("../data/BED/test.bed", "BED");
    uWriter writer("junktest.txt","GFF");

    // uWriter writer(outputName, outputType);
        while (!Parser.eof())
        {
        writer.writeToken(Parser.getNextEntry());
        }

 try
    {
        uParser GFFParser("junktest.txt", "GFF");
        while (!GFFParser.eof())
        {
            GFFParser.getNextEntry();
        }
    }
    catch(std::runtime_error& e)
    {
        cout << e.what() << endl;
    }
    catch(ugene_exception_base& e)
    {
        cout << fetchStringError(e) << endl;
    }



}
