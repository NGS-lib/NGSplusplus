#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "NGS++.h"
#include "tests_Writer_Fixtures.h"

using namespace std;
using namespace NGS;


ostringstream* m_pOsGTF = new ostringstream(ostringstream::out);
std::string gtfExpect ="chr1\t.\t.\t100\t200\t.";


TEST(TestsGTFWriter_WriteToken,TAG) {

	uTags testTag("chr1", 100, 200);
	uWriter gtfWriter(m_pOsGTF,"GTF");
    testTag.writeToOutput(gtfWriter);
	ASSERT_TRUE(m_pOsGTF->str().find(gtfExpect) != string::npos);
}

TEST(TestsGTFWriter_WriteToken,REGION) {

	uRegion testRegion("chr1", 100, 200);
	uWriter gtfWriter(m_pOsGTF,"GTF");
    testRegion.writeToOutput(gtfWriter);
	ASSERT_TRUE(m_pOsGTF->str().find(gtfExpect) != string::npos);
}

TEST(TestsGTFWriter_WriteToken,BASICNGS) {

	uBasicNGS testBasic("chr1", 100, 200);
	uWriter gtfWriter(m_pOsGTF,"GTF");
	 testBasic.writeToOutput(gtfWriter);
	ASSERT_TRUE(m_pOsGTF->str().find(gtfExpect) != string::npos);
}

TEST(TestsGTFWriter_WriteToken,GENE) {

	uGene testUgene("chr1", 100, 200);
	uWriter gtfWriter(m_pOsGTF,"GTF");
	testUgene.writeToOutput(gtfWriter);

	ASSERT_TRUE(m_pOsGTF->str().find(gtfExpect) != string::npos);
}
