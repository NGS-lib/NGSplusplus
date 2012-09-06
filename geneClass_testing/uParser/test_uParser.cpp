#include <gtest/gtest.h>
#include <string>
#include <sstream>
//TODO
/**< To test private member directly, changed private for public */
//#define private public
#include "uParser.h"

using namespace std;

/*
 * Constructor tests, filename:
 *		uParser(std::string filename, file_type type = BED);
 * 	Valid cases:
 *		ValidFileName
 *	Invalid cases:
 *		InvalidFileName
 */
TEST(uParserConstructorFilename, ValidFileName) {
	ASSERT_NO_THROW(uParser Parser("test.bed", file_type::BED));
}

TEST(uParserConstructorFilename, InvalidFileName) {
	ASSERT_THROW(uParser Parser("test2.bed", file_type::BED), std::runtime_error);
	try {
		uParser Parser("test2.bed", file_type::BED);
		ASSERT_TRUE(false);
	}
	catch (std::runtime_error& e) {
		ASSERT_STREQ(e.what(), "Error opening file: test2.bed");
	}
}

/*
 * Constructor tests, filename:
 *		uParser(std::string filename, file_type type = BED);
 * 	Valid cases:
 *		ValidStream
 *		EmptyStream
 */

TEST(uParserConstructorStream, ValidStream) {
	stringstream ss;
	ss << "chr1\t21\t31\ttest001\t.\t+\n";
	ss << "chr2\t1221\t1231\ttest002\t.\t+\n";
	ss << "chr3\t34521\t34531\ttest003\t.\t-\n";
	ss << "chr4\t456721\t456731\ttest004\t.\t-\n";
	ASSERT_NO_THROW(uParser Parser(&ss, file_type::BED));
}

TEST(uParserConstructorStream, EmptyStream) {
	stringstream ss;
	ASSERT_NO_THROW(uParser Parser(&ss, file_type::BED));
}

/*
 * Custom Constructor tests, filename:
 *		uParser(std::string filename, file_type type = BED);
 * 	Valid cases:
 *		ValidFilename
 *		AlternateDelimiter
 *	Invalid cases:
 *		InvalidFilename
 *		EmptyFieldsList
 * 		MissingMandatoryFields
 * 		CannotInferEND_POS
 */

//class uParser_test : public uParser {
//public:
//	uParser_test(const string& filename, const vector<string>& fieldsNames, char delimiter = '\t')
//		: uParser(filename, fieldsNames, delimiter) { }
//	uParser_test(istream* stream, const vector<string>& fieldsNames, char delimiter = '\t')
//		: uParser(stream, fieldsNames, delimiter) { }
//	vector<string> getFields() { return m_customFieldNames; }
//};
class CustomConstructorTests_ValidList: public ::testing::Test {
public:
	CustomConstructorTests_ValidList() {
		m_fieldsList.push_back("CHR");
		m_fieldsList.push_back("STRAND");
		m_fieldsList.push_back("START_POS");
		m_fieldsList.push_back("END_POS");
		m_fieldsList.push_back("NAME");
		m_fieldsList.push_back("JUNK");
	}
	vector<string> m_fieldsList;
};

class CustomConstructorTests_EmptyList: public ::testing::Test {
public:
	CustomConstructorTests_EmptyList() { }
	vector<string> m_fieldsList;
};

class CustomConstructorTests_InvalidListMandatoryCHR: public ::testing::Test {
public:
	CustomConstructorTests_InvalidListMandatoryCHR() { 
		m_fieldsList.push_back("STRAND");
		m_fieldsList.push_back("START_POS");
		m_fieldsList.push_back("END_POS");
		m_fieldsList.push_back("NAME");
		m_fieldsList.push_back("JUNK");
	}
	vector<string> m_fieldsList;
};

class CustomConstructorTests_InvalidListMandatorySTART_POS: public ::testing::Test {
public:
	CustomConstructorTests_InvalidListMandatorySTART_POS() { 
		m_fieldsList.push_back("CHR");
		m_fieldsList.push_back("STRAND");
		m_fieldsList.push_back("END_POS");
		m_fieldsList.push_back("NAME");
		m_fieldsList.push_back("JUNK");
	}
	vector<string> m_fieldsList;
};

class CustomConstructorTests_InvalidListEND_POS: public ::testing::Test {
public:
	CustomConstructorTests_InvalidListEND_POS() {
		m_fieldsList.push_back("CHR");
		m_fieldsList.push_back("STRAND");
		m_fieldsList.push_back("START_POS");
		m_fieldsList.push_back("NAME");
		m_fieldsList.push_back("JUNK");
	}
	vector<string> m_fieldsList;
};

TEST_F(CustomConstructorTests_ValidList, Filename_ValidFilename) {
	ASSERT_NO_THROW(uParser Parser("test.custom", m_fieldsList));
}

TEST_F(CustomConstructorTests_ValidList, Filename_AlternateDelimiter) {
	ASSERT_NO_THROW(uParser Parser("test.csv", m_fieldsList, ','));
}

TEST_F(CustomConstructorTests_ValidList, Filename_InvalidFilename) {
	ASSERT_THROW(uParser Parser("asdf", m_fieldsList), runtime_error);
}

TEST_F(CustomConstructorTests_EmptyList, ValidFilename) {
	ASSERT_THROW(uParser Parser("test.custom", m_fieldsList), customParser_missing_mandatory_values);
}

TEST_F(CustomConstructorTests_InvalidListMandatoryCHR, ValidFilename) {
	ASSERT_THROW(uParser Parser("test.custom", m_fieldsList), customParser_missing_mandatory_values);
}

TEST_F(CustomConstructorTests_InvalidListMandatorySTART_POS, ValidFilename) {
	ASSERT_THROW(uParser Parser("test.custom", m_fieldsList), customParser_missing_mandatory_values);
}

TEST_F(CustomConstructorTests_InvalidListEND_POS, ValidFilename) {
	ASSERT_THROW(uParser Parser("test.custom", m_fieldsList), customParser_missing_mandatory_values);
}

/*
 * Custom Constructor tests, stream:
 *		uParser(std::string filename, file_type type = BED);
 * 	Valid cases:
 *		ValidStream
 *		AlternateDelimiter
 *		EmptyStream	
 *	Invalid cases:
 *		EmptyFieldsList
 * 		MissingMandatoryFields
 * 		CannotInferEND_POS
 */

TEST_F(CustomConstructorTests_ValidList, Stream_ValidStream) {
	stringstream ss;
	ss << "chr1\t21\t31\ttest001\t.\t+\n";
	ss << "chr2\t1221\t1231\ttest002\t.\t+\n";
	ss << "chr3\t34521\t34531\ttest003\t.\t-\n";
	ss << "chr4\t456721\t456731\ttest004\t.\t-\n";
	ASSERT_NO_THROW(uParser Parser(&ss, m_fieldsList));
}

TEST_F(CustomConstructorTests_ValidList, Stream_EmptyStream) {
	stringstream ss;
	ASSERT_NO_THROW(uParser Parser(&ss, m_fieldsList));
}

TEST_F(CustomConstructorTests_ValidList, Stream_AlternateDelimiter) {
	stringstream ss;
	ASSERT_NO_THROW(uParser Parser(&ss, m_fieldsList, ','));
}

TEST_F(CustomConstructorTests_EmptyList, ValidStream) {
	stringstream ss;
	ASSERT_THROW(uParser Parser(&ss, m_fieldsList), customParser_missing_mandatory_values);
}

TEST_F(CustomConstructorTests_InvalidListMandatoryCHR, ValidStream) {
	stringstream ss;
	ASSERT_THROW(uParser Parser(&ss, m_fieldsList), customParser_missing_mandatory_values);
}

TEST_F(CustomConstructorTests_InvalidListMandatorySTART_POS, ValidStream) {
	stringstream ss;
	ASSERT_THROW(uParser Parser(&ss, m_fieldsList), customParser_missing_mandatory_values);
}

TEST_F(CustomConstructorTests_InvalidListEND_POS, ValidStream) {
	stringstream ss;
	ASSERT_THROW(uParser Parser(&ss, m_fieldsList), customParser_missing_mandatory_values);
}

/* Tests for the function:
 * 		bool eof() const 
 *	Valid cases:
 *		NotEndOfFile
 *		EndOfFile
 *		NoTokenInStream 
 */

TEST(uParserEof, NotEndOfFile) {
	uParser Parser("test.bed", file_type::BED);
	ASSERT_FALSE(Parser.eof());
	uToken Token = Parser.getNextEntry();
	ASSERT_FALSE(Parser.eof());
}

TEST(uParserEof, EndOfFile) {
	uParser Parser("test.bed", file_type::BED);
	uToken Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	ASSERT_TRUE(Parser.eof());
}

TEST(uParserEof, NoTokenInStream) {
	stringstream ss;
	uParser Parser(&ss, file_type::BED);
	ASSERT_TRUE(Parser.eof());
}

/* 
 * Tests for the function:
 *		uToken getNextEntry();
 *	Valid cases:
 *		CorrectlyFormatedBED6
 *		CorrectlyFormatedBED4
 *		CorrectlyFormatedCustom
 *		CorrectlyFormatedCustomAlternateDelimiter
 *		
 *	Invalid cases:
 *		IncorrectlyFormatedBED
 *		IncorrectlyFormatedCustom
 *		ReachedEOF
 */

TEST(uParserGetNextEntry, CorrectlyFormatedBED6) {
	uParser Parser("test.bed", file_type::BED);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "test001");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");

	Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1221");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "1231");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "test002");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}

TEST(uParserGetNextEntry, CorrectlyFormatedBED4) {
	uParser Parser("test.bed", file_type::BED);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "test001");

	Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1221");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "1231");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "test002");
}

TEST_F(CustomConstructorTests_ValidList, getNextEntry_CorrectlyFormatedCustom) {
	uParser Parser("test.custom", m_fieldsList);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");

	Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1221");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "1231");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}

TEST_F(CustomConstructorTests_ValidList, getNextEntry_CorrectlyFormatedCustomAlternateDelimiter) {
	uParser Parser("test.csv", m_fieldsList, ',');
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");

	Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1221");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "1231");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}

TEST(uParserGetNextEntry, IncorrectlyFormatedBED) {
	uParser Parser("incorrect.bed", file_type::BED);
	ASSERT_THROW(Parser.getNextEntry(), invalid_value_throw);

	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1221");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "1231");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "test002");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}

TEST_F(CustomConstructorTests_ValidList, getNextEntry_IncorrectlyFormatedCustom) {
	uParser Parser("incorrect.custom", m_fieldsList);
	ASSERT_THROW(Parser.getNextEntry(), invalid_value_throw);
	ASSERT_THROW(Parser.getNextEntry(), invalid_uToken_throw);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr3");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "34521");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "34531");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "-");
}

TEST(uParserGetNextEntry, ReachedEOF) {
	uParser Parser("test.bed", file_type::BED);
	uToken Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	ASSERT_THROW(Token = Parser.getNextEntry(), end_of_file_throw);
}
