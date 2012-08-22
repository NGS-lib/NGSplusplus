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
