#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "NGS++.h"

using namespace std;
using namespace NGS;


const std::string expectedString="chr1	10301	10321	4.0\nchr1	10321	10326	4.1\nchr1	10326	10331	4.2\nchr1	10331	10336	4.3\nchr1	10336	10345	4.4\nchr1	10345	10348	4.5\nchr1	10348	10357	4.4\nchr1	10357	10362	4.3\nchr1	10362	10367	4.2\nchr1	10367	10371	4.1\n";

TEST(TestsBedGraphWriter, ReadWriteLines) {
    uParser ourParser("../data/bedgraph/Small10K.bedgraph", "BEDGRAPH");
    auto m_pOssBedGraph = new ostringstream(ostringstream::out);
    uWriter ourWriter(m_pOssBedGraph, "BEDGRAPH");
    for (int i=0; i<10; i++)
        ourWriter.writeToken(ourParser.getNextEntry());

    EXPECT_EQ(m_pOssBedGraph->str(),expectedString);
}



ostringstream* m_pOsBedGraph = new ostringstream(ostringstream::out);
std::string bedGraph ="chr1\t100\t200\t0\n";

TEST(TestsbedGraphWriter_WriteToken,TAG) {

	uTags testTag("chr1", 100, 200);
	uWriter bedGraphWriter(m_pOsBedGraph,"BEDGRAPH");
    testTag.writeToOutput(bedGraphWriter);
	ASSERT_TRUE(m_pOsBedGraph->str().find(bedGraph) != string::npos);
}

TEST(TestsbedGraphWriter_WriteToken,REGION) {

	uRegion testRegion("chr1", 100, 200);
	uWriter bedGraphWriter(m_pOsBedGraph,"BED");
    testRegion.writeToOutput(bedGraphWriter);
	ASSERT_TRUE(m_pOsBedGraph->str().find(bedGraph) != string::npos);
}

TEST(TestsbedGraphWriter_WriteToken,BASICNGS) {

	uBasicNGS testBasic("chr1", 100, 200);
	uWriter bedGraphWriter(m_pOsBedGraph,"BED");
	 testBasic.writeToOutput(bedGraphWriter);
	ASSERT_TRUE(m_pOsBedGraph->str().find(bedGraph) != string::npos);
}

TEST(TestsbedGraphWriter_WriteToken,GENE) {

	uGene testUgene("chr1", 100, 200);
	uWriter bedGraphWriter(m_pOsBedGraph,"BED");
	testUgene.writeToOutput(bedGraphWriter);

	ASSERT_TRUE(m_pOsBedGraph->str().find(bedGraph) != string::npos);
}

