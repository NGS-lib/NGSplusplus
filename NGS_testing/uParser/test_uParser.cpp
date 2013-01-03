#include <gtest/gtest.h>
#include <string>
#include <sstream>
//TODO : Split tests by type of file??
/**< To test private member directly, changed private for public */
//#define private public
#include "NGS++.h"
#include "fixtures.h"

using namespace std;
using namespace NGS;

/* Tests for the function:
 * 		bool eof() const
 *	Valid cases:
 *		NotEndOfFile
 *		EndOfFile
 *		NoTokenInStream
 */
TEST(uParserEof, NotEndOfFile) {
	uParser Parser("../data/BED/test.bed", "BED");
	ASSERT_FALSE(Parser.eof());
	uToken Token = Parser.getNextEntry();
	ASSERT_FALSE(Parser.eof());
}

TEST(uParserEof, EndOfFile) {
	uParser Parser("../data/BED/test.bed", "BED");
	uToken Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	Token = Parser.getNextEntry();
	ASSERT_TRUE(Parser.eof());
}

TEST(uParserEof, NoTokenInStream) {
	stringstream ss;
	uParser Parser(&ss, "BED");
	ASSERT_TRUE(Parser.eof());
}

/*
 * Tests for the function:
 *		std::string getUnformatedHeader() const;
 *	Valid cases:
 *		NoHeader
 *		WithHeader
 */

TEST(getUnformatedHeader, NoHeader) {
	uParser Parser("../data/BED/test.bed", "BED");
	ASSERT_EQ(Parser.getUnformatedHeader(), "");
}

TEST(getUnformatedHeader, WithHeader) {
	uParser Parser("../data/BED/header.bed", "BED", true);
	string unformated = "";
	unformated += "browser position chr7:127471196-127495720\n";
	unformated += "browser hide all\n";
	unformated += "track name=\"ItemRGBDemo\" description=\"Item RGB demonstration\" visibility=2\n";
	unformated += "itemRgb=\"On\"\n";
	ASSERT_EQ(Parser.getUnformatedHeader(), unformated);
}
