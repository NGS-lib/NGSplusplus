#include <gtest/gtest.h>
#include "uParser.h"
#include <string>
#include <sstream>

using namespace std;

/*
 * Constructor tests:
 *		uParser(std::string filename, file_type type = BED);
 * 	Valid cases:
 *		ValidFileName
 *	Invalid cases:
 *		InvalidFileName
 */
TEST(uParserConstructor, ValidFileName) {
	ASSERT_NO_THROW(uParser Parser("test.bed", file_type::BED));
}

TEST(uParserConstructor, InvalidFileName) {
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
 * Test for the function:
 *		uToken getNextEntry();
 *	Valid cases:
 *		CorrectlyFormatedBED6
 *		CorrectlyFormatedBED4
 *		
 *	Invalid cases:
 *		IncorrectlyFormatedBED
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

TEST(uParserGetNextEntry, IncorrectlyFormatedBED) {
	uParser Parser("incorrect.bed", file_type::BED);
	ASSERT_THROW(Parser.getNextEntry(), invalid_uToken_throw);

	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1221");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "1231");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "test002");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}

TEST(uParserGetNextEntry, ReachedEOF) {
	uParser Parser("test.bed", file_type::BED);
	uToken Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	ASSERT_THROW(Token = Parser.getNextEntry(), end_of_file_throw);
}
