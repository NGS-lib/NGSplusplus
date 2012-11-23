#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <sstream>
#include "NGS++.h"

using namespace std;
using namespace NGS;

class validTokens {
public:
	validTokens() {
		/**< Token 1: Contain all token_param */
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
        /**< Token 5: Simple token, START/END CHR only */
	//	ss.clear();
	//	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	//	uToken Token6(ss);
	//	m_vTokens.push_back(Token6);
	}

	vector<uToken> m_vTokens;
};



class TestGFFWriter: public ::testing::Test {
public:
	TestGFFWriter() {
		/**< Initialize streams */
        m_pOssSAM = new ostringstream(ostringstream::out);
		/**< Initialize writers */
		uWriter writerGFF(m_pOssSAM, "GFF");
		/**< Write tokens */
		validTokens tokens;
		for (size_t i = 0; i < tokens.m_vTokens.size(); i++) {
			writerGFF.writeToken(tokens.m_vTokens[i]);
		}
	}
	~TestGFFWriter() {
		delete m_pOssSAM;
		m_pOssSAM = nullptr;

	}
	ostringstream* m_pOssSAM = nullptr;
};




class TestSamWriter: public ::testing::Test {
public:
	TestSamWriter() {
		/**< Initialize streams */
        m_pOssSAM = new ostringstream(ostringstream::out);
		/**< Initialize writers */
		uWriter writerSAM(m_pOssSAM, "SAM");

		/**< Write tokens */
		validTokens tokens;
		for (size_t i = 0; i < tokens.m_vTokens.size(); i++) {
			writerSAM.writeToken(tokens.m_vTokens[i]);
		}
	}
	~TestSamWriter() {
		delete m_pOssSAM;
		m_pOssSAM = nullptr;

	}
	ostringstream* m_pOssSAM = nullptr;
};


class TestsBedWriter: public ::testing::Test {
public:
	TestsBedWriter() {
		/**< Initialize streams */
		m_pOssBed4 = new ostringstream(ostringstream::out);
		m_pOssBed6 = new ostringstream(ostringstream::out);
		/**< Initialize writers */
		uWriter writerBed4(m_pOssBed4, "BED4");
		uWriter writerBed6(m_pOssBed6, "BED6");
		/**< Write tokens */
		validTokens tokens;
		for (size_t i = 0; i < tokens.m_vTokens.size(); i++) {
			writerBed4.writeToken(tokens.m_vTokens[i]);
			writerBed6.writeToken(tokens.m_vTokens[i]);
		}
	}

	~TestsBedWriter() {
		delete m_pOssBed4;
		m_pOssBed4 = nullptr;
		delete m_pOssBed6;
		m_pOssBed6 = nullptr;
	}

	ostringstream* m_pOssBed4 = nullptr;
	ostringstream* m_pOssBed6 = nullptr;
};

class TestsCustomWriter: public ::testing::Test {
public:
	TestsCustomWriter() {
		/**< Initialize fields list */
		m_fieldList.push_back("CHR");
		m_fieldList.push_back("START_POS");
		m_fieldList.push_back("SEQUENCE");
		/**< Initialize streams */
		m_pOss = new ostringstream(ostringstream::out);
		/**< Initialize writers */
		uWriter writerCustom(m_pOss, m_fieldList);
		/**< Write tokens */
		validTokens tokens;
		for (size_t i = 0; i < tokens.m_vTokens.size(); i++) {
			writerCustom.writeToken(tokens.m_vTokens[i]);
		}
	}

	vector<string> m_fieldList;
	ostringstream* m_pOss = nullptr;
};

class TestsCustomWriterOneInvalidValue: public ::testing::Test {
public:
	TestsCustomWriterOneInvalidValue() {
		/**< Initialize fields list */
		m_fieldList.push_back("CHR");
		m_fieldList.push_back("START_POS");
		m_fieldList.push_back("SEQUENCE");
		m_fieldList.push_back("NEW_SCORE");
		/**< Initialize streams */
		m_pOss = new ostringstream(ostringstream::out);
		/**< Initialize writers */
		uWriter writerCustom(m_pOss, m_fieldList);
		/**< Write tokens */
		validTokens tokens;
		for (size_t i = 0; i < tokens.m_vTokens.size(); i++) {
			writerCustom.writeToken(tokens.m_vTokens[i]);
		}
	}

	vector<string> m_fieldList;
	ostringstream* m_pOss = nullptr;
};

class TestsCustomWriterNoValidValues: public ::testing::Test {
public:
	TestsCustomWriterNoValidValues() {
		/**< Initialize fields list */
		m_fieldList.push_back("NEW_SCORE");
		m_fieldList.push_back("OLD_SCORE");
		m_fieldList.push_back("DENSITY");
		/**< Initialize streams */
		m_pOss = new ostringstream(ostringstream::out);
		/**< Initialize writers */
		uWriter writerCustom(m_pOss, m_fieldList);
		/**< Write tokens */
		validTokens tokens;
		for (size_t i = 0; i < tokens.m_vTokens.size(); i++) {
			writerCustom.writeToken(tokens.m_vTokens[i]);
		}
	}

	vector<string> m_fieldList;
	ostringstream* m_pOss = nullptr;
};
