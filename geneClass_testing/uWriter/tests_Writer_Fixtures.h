#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "IO/Writer/uWriter.h"

using namespace std;
using namespace NGS;

class TestsBedWriter: public ::testing::Test {
public:
	TestsBedWriter() {
		// Initialize write
		m_pOss = new ostringstream(ostringstream::out);	
		m_pWriter = new uWriter(m_pOss, "BED");
		// Write a token with writer
		stringstream ss;
		ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
		ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
		ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
		ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
		uToken Token1(ss);
		m_pWriter->writeToken(Token1);
		// Write a second token with writer
		ss.clear();
		ss << "CHR\tchr2\n" << "START_POS\t101\n" << "END_POS\t121\n";
		ss << "STRAND\t-\n" << "MAP_SCORE\t24\n" << "PHRED_SCORE\t####################\n";
		ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tNCGTN.acgtn.ACGTGTCN\n";
		ss << "SEQ_NAME\tab00002\n" << "FLAGS\t257\n";
		uToken Token2(ss);
		m_pWriter->writeToken(Token2);

	}
	~TestsBedWriter() {
		delete m_pOss;
		m_pOss = nullptr;
		delete m_pWriter;
		m_pWriter = nullptr;
	}
	ostringstream* m_pOss = nullptr;
	uWriter* m_pWriter = nullptr;
};
