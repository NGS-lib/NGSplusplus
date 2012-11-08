#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "IO/Parser/uParser.h"
//#include "fixtures.h"

using namespace std;
using namespace NGS;

/*
 * Constructor tests, filename:
 *		uParser(const std::string& filename, const std::string & type, bool header = false);
 * 	Valid cases:
 *		ValidFileName
 *		CorrectlyFormatedHeaderBed
 *	Invalid cases:
 *		InvalidFileName
 *		IncorrectlyFormatedHeader //TODO
 */
TEST(uParserConstructorFilename, ValidFileName) {

	ASSERT_NO_THROW(uParser Parser("../data/BED/test.bed", "BED"));
}

TEST(uParserConstructorFilename, CorrectlyFormatedHeaderBed) {
	ASSERT_NO_THROW(uParser Parser("../data/BED/header.bed", "BED", true));
}

TEST(uParserConstructorFilename, InvalidFileName) {
	ASSERT_THROW(uParser Parser("../data/BED/test2.bed", "BED"), std::runtime_error);
	try {
		uParser Parser("../data/BED/test2.bed", "BED");
		ASSERT_TRUE(false);
	}
	catch (std::runtime_error& e) {
		ASSERT_STREQ(e.what(), "Error opening file: ../data/BED/test2.bed");
	}
}

/*
 * Constructor tests, stream:
 *		uParser(std::iostream* stream, const std::string & type, bool header = false);
 * 	Valid cases:
 *		ValidStream
 *		CorrectlyFormatedHeader
 *		EmptyStream
 *	Invalid case:
 *		IncorrectlyFormatedHeader //TODO
 */

TEST(uParserConstructorStream, ValidStream) {
	stringstream ss;
	ss << "chr1\t21\t31\ttest001\t.\t+\n";
	ss << "chr2\t1221\t1231\ttest002\t.\t+\n";
	ASSERT_NO_THROW(uParser Parser(&ss, "BED"));
}

TEST(uParserConstructorStream, CorrectlyFormatedHeader) {
	stringstream ss;
	ss << "browser position chr7:127471196-127495720\n";
	ss << "chr1\t21\t31\ttest001\t.\t+\n";
	ss << "chr2\t1221\t1231\ttest002\t.\t+\n";
	ASSERT_NO_THROW(uParser Parser(&ss, "BED", true));
}

TEST(uParserConstructorStream, EmptyStream) {
	stringstream ss;
	ASSERT_NO_THROW(uParser Parser(&ss, "BED"));
}

/*
 * Tests for the function:
 *		uToken getNextEntry();
 *	Valid cases:
 *		CorrectlyFormatedBED6
 *		CorrectlyFormatedBED4
 *		CorrectlyFormatedHeaderBED
 *
 *	Invalid cases:
 *		IncorrectlyFormatedBED
 * 		IncorrectlyFormatedHeaderBED
 *		CorrectlyFormatedHeaderButNotSpecifiedBED
 *		ReachedEOF
 */
// TODO: Check score also!
TEST(uParserGetNextEntry, CorrectlyFormatedBED6) {
	uParser Parser("../data/BED/test.bed", "BED");
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
	uParser Parser("../data/BED/test.bed", "BED");
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

TEST(uParserGetNextEntry, CorrectlyFormatedHeaderBED) {
	uParser Parser("../data/BED/header.bed", "BED", true);
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

TEST(uParserGetNextEntry, IncorrectlyFormatedBED) {
	uParser Parser("../data/BED/incorrect.bed", "BED", false);
	ASSERT_THROW(Parser.getNextEntry(), invalid_value_throw);

	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1221");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "1231");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "test002");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}
// TODO: Check if next entry is ok
TEST(uParserGetNextEntry, IncorrectlyFormatedHeaderBED) {
	uParser Parser("../data/BED/incorrect_header.bed", "BED", true);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
	ASSERT_THROW(Parser.getNextEntry(), invalid_value_throw);
}

TEST(uParserGetNextEntry, CorrectlyFormatedHeaderButNotSpecifiedBED) {
	uParser Parser("../data/BED/header.bed", "BED");
	ASSERT_THROW(Parser.getNextEntry(), ugene_exception_base);
	ASSERT_THROW(Parser.getNextEntry(), ugene_exception_base);
	ASSERT_THROW(Parser.getNextEntry(), ugene_exception_base);
	ASSERT_THROW(Parser.getNextEntry(), ugene_exception_base);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}

TEST(uParserGetNextEntry, ReachedEOF) {
	uParser Parser("../data/BED/test.bed", "BED");
	uToken Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	ASSERT_THROW(Token = Parser.getNextEntry(), end_of_file_throw);
}
