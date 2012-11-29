#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "NGS++.h"
#include "tests_Writer_Fixtures.h"

using namespace std;
using namespace NGS;


TEST(TestsWigWrite, WriteValid) {
    /**< Initialize streams */
   //  std::ostream & outFile = ((outputPath.size()!=0) ? outputOS : std::cout);
        vector<uToken> m_vTokens;
		stringstream ss;
		ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
		ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
		ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
		ss << "SEQ_NAME\tab00001\n" << "SCORE\t111\n" << "FLAGS\t256\n";
		uToken Token1(ss);
		m_vTokens.push_back(Token1);
		/**< Token 2: Token 1 minus SCORE */
		ss.clear();
		ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
		ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
		ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
		ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
		uToken Token2(ss);
		m_vTokens.push_back(Token2);
		/**< Token 3: Token 1 minus SEQ_NAME */
		ss.clear();
		ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
		ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
		ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
		ss << "SCORE\t111\n" << "FLAGS\t256\n";
		uToken Token3(ss);
		m_vTokens.push_back(Token3);
		/**< Token 4: Token 1 minus SCORE and STRAND */
		ss.clear();
		ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
		ss << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
		ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
		ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
		uToken Token4(ss);
		m_vTokens.push_back(Token4);
		/**< Token 5: Token 1 minus SEQUENCE */
		ss.clear();
		ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
		ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
		ss << "CIGAR\t2M3I2X1=12X\n";
		ss << "SEQ_NAME\tab00001\n" << "SCORE\t111\n" << "FLAGS\t256\n";
		uToken Token5(ss);
		m_vTokens.push_back(Token5);

     //   std::ostream& coutStream=std::cout;

}

TEST(TestsWigWrite, ParseBedWriteWig)
{
    std::ostream& coutStream=std::cout;
    uParser Parser("../data/BED/header.bed", "BED", true);
    uWriter wigWrite(&coutStream,"WIG");
    while(!(Parser.eof()))
        ASSERT_NO_THROW(wigWrite.writeToken(Parser.getNextEntry()));

}

TEST(TestsWigWrite, ParseSamWriteWig)
{
    std::ostream& coutStream=std::cout;
    uParser ourParser("../data/SAM/fiveCountValid.sam", "SAM");
    uWriter wigWrite(&coutStream,"WIG");
    while(!(ourParser.eof()))
        wigWrite.writeToken(ourParser.getNextEntry());

	//ASSERT_NO_THROW(uParser Parser("../data/BED/header.bed", "BED", true));
}



