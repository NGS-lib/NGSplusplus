#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "NGS++.h"
//#include "fixtures.h"

using namespace std;
using namespace NGS;


const std::string MACSPATH="../data/BED/macs_test.bed";
const std::string ENDLPATH="../data/BED/extrEndl.bed";
const std::string INCORRECT_HEADER="../data/BED/incorrect_header.bed";
const std::string HEADER="../data/BED/header.bed";
const std::string INCORRECT="../data/BED/incorrect.bed";
const std::string TEST="../data/BED/test.bed";
const std::string BED3="../data/BED/bed3.bed";
const std::string BED5="../data/BED/bed5.bed";


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
TEST(uParserBedConstructorFilename, ValidFileName) {

	ASSERT_NO_THROW(uParser Parser(TEST, "BED"));
}

TEST(uParserBedGetNextEntry, ExtraENDL) {

    uParser Parser(ENDLPATH, "BED");
	uToken Token = Parser.getNextEntry();
	ASSERT_NO_THROW(while(!(Parser.eof()))
                    Parser.getNextEntry();
                 );
}

TEST(uParserBedConstructorFilename, CorrectlyFormatedHeaderBed) {
	ASSERT_NO_THROW(uParser Parser(HEADER, "BED", true));
}

TEST(uParserBedConstructorFilename, InvalidFileName) {
	ASSERT_THROW(uParser Parser("../data/BED/notHere.bed", "BED"), std::runtime_error);
	try {
		uParser Parser("../data/BED/notHere.bed", "BED");
		ASSERT_TRUE(false);
	}
	catch (std::runtime_error& e) {
		ASSERT_STREQ(e.what(), "Error opening file: ../data/BED/notHere.bed");
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

TEST(uParserBedConstructorStream, ValidStream) {
	stringstream ss;
	ss << "chr1\t21\t31\ttest001\t.\t+\n";
	ss << "chr2\t1221\t1231\ttest002\t.\t+\n";
	ASSERT_NO_THROW(uParser Parser(&ss, "BED"));
}

TEST(uParserBedConstructorStream, CorrectlyFormatedHeader) {
	stringstream ss;
	ss << "browser position chr7:127471196-127495720\n";
	ss << "chr1\t21\t31\ttest001\t.\t+\n";
	ss << "chr2\t1221\t1231\ttest002\t.\t+\n";
	ASSERT_NO_THROW(uParser Parser(&ss, "BED", true));
}

TEST(uParserBedConstructorStream, EmptyStream) {
	stringstream ss;
	ASSERT_NO_THROW(uParser Parser(&ss, "BED"));
}


/**< BlackBox data tests from applications */


 TEST(uParserBed_GetNextEntry, MACSREAL) {
	uParser Parser(MACSPATH, "BED");
    while(Parser.eof()==false)
    {
        EXPECT_NO_THROW(Parser.getNextEntry());
    }
}

 TEST(uParserBed_GetNextEntry, BED3) {
	uParser Parser(BED3, "BED");
    while(Parser.eof()==false)
    {
        EXPECT_NO_THROW(Parser.getNextEntry());
    }
}

 TEST(uParserBed_GetNextEntry, BED5) {
	uParser Parser(BED5, "BED");
    while(Parser.eof()==false)
    {
        EXPECT_NO_THROW(Parser.getNextEntry());
    }
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
TEST(uParserBedGetNextEntry, CorrectlyFormatedBED6) {
	uParser Parser(TEST, "BED");
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

TEST(uParserBedGetNextEntry, CorrectlyFormatedBED4) {
	uParser Parser(TEST, "BED");
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

TEST(uParserBedGetNextEntry, CorrectlyFormatedHeaderBED) {
	uParser Parser(HEADER, "BED", true);
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

TEST(uParserBedGetNextEntry, IncorrectlyFormatedBED) {
	uParser Parser(INCORRECT, "BED", false);
	ASSERT_THROW(Parser.getNextEntry(), invalid_value_throw);

	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1221");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "1231");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "test002");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
}
// TODO: Check if next entry is ok
TEST(uParserBedGetNextEntry, IncorrectlyFormatedHeaderBED) {
	uParser Parser(INCORRECT_HEADER, "BED", true);
	uToken Token = Parser.getNextEntry();
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "31");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
	ASSERT_THROW(Parser.getNextEntry(), invalid_value_throw);
}

TEST(uParserBedGetNextEntry, CorrectlyFormatedHeaderButNotSpecifiedBED) {
	uParser Parser(HEADER, "BED");
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

TEST(uParserBedGetNextEntry, ReachedEOF) {
	uParser Parser(TEST, "BED");
	uToken Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	ASSERT_THROW(Token = Parser.getNextEntry(), end_of_file_throw);
}


TEST(uParserBedGetNextEntry, MACSSAMPLE) {
	uParser Parser(HEADER, "BED", true);
}

