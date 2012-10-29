#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "IO/Parser/uParser.h" 
#include "fixtures.h"

using namespace std;
using namespace NGS;

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

TEST_F(CustomConstructorTests_ValidList, Filename_ValidFilename) {
	ASSERT_NO_THROW(uParser Parser("test.custom", m_fieldsList));
}

TEST_F(CustomConstructorTests_ValidList, Filename_AlternateDelimiter) {
	ASSERT_NO_THROW(uParser Parser("test.csv", m_fieldsList, false, ','));
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
	ASSERT_NO_THROW(uParser Parser(&ss, m_fieldsList, false, ','));
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

/*
 * Tests for the function:
 *		uToken getNextEntry();
 *	Valid cases:
 *		CorrectlyFormatedHeaderCUSTOM
 *		CorrectlyFormatedCustomAlternateDelimiter
 *
 *	Invalid cases:
 *		IncorrectlyFormatedCustom
 * 		IncorrectlyFormatedHeaderCUSTOM
 *		CorrectlyFormatedHeaderButNotSpecifiedCUSTOM
 */

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

TEST_F(CustomConstructorTests_ValidList, getNextEntry_CorrectlyFormatedHeaderCustom) {
	uParser Parser("header.custom", m_fieldsList, true);
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
	uParser Parser("test.csv", m_fieldsList, false, ',');
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

TEST_F(CustomConstructorTests_ValidList, getNextEntry_IncorrectlyFormatedCustom) {
	uParser Parser("incorrect.custom", m_fieldsList, false);
	ASSERT_THROW(Parser.getNextEntry(), invalid_value_throw);
	ASSERT_THROW(Parser.getNextEntry(), invalid_uToken_throw);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr3");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "34521");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "34531");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "-");
}

TEST_F(CustomConstructorTests_ValidList, getNextEntry_IncorrectlyFormatedHeaderCUSTOM) {
	uParser Parser("incorrect_header.custom", m_fieldsList, true);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
	ASSERT_THROW(Parser.getNextEntry(), uToken_exception_base);
}

TEST_F(CustomConstructorTests_ValidList, getNextEntry_CorrectlyFormatedHeaderButNotSpecifiedBED) {
	uParser Parser("header.custom", m_fieldsList);
	ASSERT_THROW(Parser.getNextEntry(), uToken_exception_base);
	ASSERT_THROW(Parser.getNextEntry(), uToken_exception_base);
	ASSERT_THROW(Parser.getNextEntry(), uToken_exception_base);
	ASSERT_THROW(Parser.getNextEntry(), uToken_exception_base);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}
