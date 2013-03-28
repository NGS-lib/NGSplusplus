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
	cout << m_pOssGFF->str();
    ASSERT_TRUE(m_pOssGFF->str().find(expected) == string::npos);
}
TEST_F(TestGFFWriter, WriteGFFCol3) {
	string expected = "chr1	.	.	1	21	111	+	.";
    ASSERT_TRUE(m_pOssGFF->str().find(expected) != string::npos);
}

TEST(TestGFFWriterConvert, FromBedReadGFF) {
   //  auto p_SSGFF = new ostringstream(ostringstream::out);
    uParser Parser("../data/BED/test.bed", "BED");
    uWriter writer("junktest.txt","UCSCGFF");

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


ostringstream* m_pOsGFF = new ostringstream(ostringstream::out);
std::string gffExpect ="chr1\t.\t.\t100\t200\t.";

TEST(TestsGFFWriter_WriteToken,TAG) {

	uTags testTag("chr1", 100, 200);
	uWriter gffWriter(m_pOsGFF,"UCSCGFF");
    testTag.writeToOutput(gffWriter);
	ASSERT_TRUE(m_pOsGFF->str().find(gffExpect) != string::npos);
}

TEST(TestsGFFWriter_WriteToken,REGION) {

	uRegion testRegion("chr1", 100, 200);
	uWriter gffWriter(m_pOsGFF,"UCSCGFF");
    testRegion.writeToOutput(gffWriter);
	ASSERT_TRUE(m_pOsGFF->str().find(gffExpect) != string::npos);
}

TEST(TestsGFFWriter_WriteToken,BASICNGS) {

	uBasicNGS testBasic("chr1", 100, 200);
	uWriter gffWriter(m_pOsGFF,"UCSCGFF");
	 testBasic.writeToOutput(gffWriter);
	ASSERT_TRUE(m_pOsGFF->str().find(gffExpect) != string::npos);
}

TEST(TestsGFFWriter_WriteToken,GENE) {

	uGene testUgene("chr1", 100, 200);
	uWriter gffWriter(m_pOsGFF,"UCSCGFF");
	testUgene.writeToOutput(gffWriter);

	ASSERT_TRUE(m_pOsGFF->str().find(gffExpect) != string::npos);
}













