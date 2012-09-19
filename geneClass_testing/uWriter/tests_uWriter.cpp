#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "IO/Writer/uWriter.h"
#include "tests_Writer_Fixtures.h"

using namespace std;
using namespace NGS;

//TODO: Ajouter les tests pour InvalidFileType quand ce sera implémenté
/*
 * Tests for uWriter filename constructor: 
 * 		uWriter(const std::string& filename, const std::string& type);
 *	Valid case:
 *		NormalUsage
 *	Invalid case:
 *		EmptyFilename
 */

TEST(uWriterFilenameConstructor, NormalUsage) {
	ASSERT_NO_THROW(uWriter writer("test.txt", "BED4"));
	ASSERT_NO_THROW(uWriter writer("test.txt", "BED6"));
	ASSERT_NO_THROW(uWriter writer("test.txt", "CUSTOM"));
}

TEST(uWriterFilenameConstructor, EmptyFilename) {
	ASSERT_THROW(uWriter writer("", "BED4"), std::runtime_error);
	ASSERT_THROW(uWriter writer("", "BED6"), std::runtime_error);
	ASSERT_THROW(uWriter writer("", "CUSTOM"), std::runtime_error);
}

/*
 * Tests for uWriter ostream constructor: 
 * 		uWriter(std::ostream* os, const std::string& type);
 *	Valid case:
 *		NormalUsage
 *	Invalid case:
 *		InvalidOstream
 */

TEST(uWriterOstreamConstructor, NormalUsage) {
	ostringstream oss;
	ASSERT_NO_THROW(uWriter writer(&oss, "BED4"));
	ASSERT_NO_THROW(uWriter writer(&oss, "BED6"));
	ASSERT_NO_THROW(uWriter writer(&oss, "CUSTOM"));
}

TEST(uWriterOstreamConstructor, InvalidOstream) {
	ostringstream* oss = nullptr;
	ASSERT_THROW(uWriter writer(oss, "BED4"), std::runtime_error);
	ASSERT_THROW(uWriter writer(oss, "BED6"), std::runtime_error);
	ASSERT_THROW(uWriter writer(oss, "CUSTOM"), std::runtime_error);
}

/*
 * Tests for uWriter filename custom constructor: 
 * 		uWriter(const std::string& filename, const std::vector<std::string>& fieldsNames, const std::string& type);
 *	Valid case:
 *		NormalUsage
 *	Invalid case:
 *		EmptyFilename
 *		EmptyFieldsNames
 */

TEST(uWriterCustomFilenameConstructor, NormalUsage) {
	vector<string> vs;
	vs.push_back("CHR");
	ASSERT_NO_THROW(uWriter writer("test.txt", vs, "BED4"));
	ASSERT_NO_THROW(uWriter writer("test.txt", vs, "BED6"));
	ASSERT_NO_THROW(uWriter writer("test.txt", vs, "CUSTOM"));
}

TEST(uWriterCustomFilenameConstructor, EmptyFilename) {
	vector<string> vs;
	vs.push_back("CHR");
	ASSERT_THROW(uWriter writer("", vs, "BED4"), std::runtime_error);
	ASSERT_THROW(uWriter writer("", vs, "BED6"), std::runtime_error);
	ASSERT_THROW(uWriter writer("", vs, "CUSTOM"), std::runtime_error);
}

TEST(uWriterCustomFilenameConstructor, EmptyFieldsNames) {
	vector<string> vs;
	ASSERT_THROW(uWriter writer("test.txt", vs, "BED4"), no_fields_names);
	ASSERT_THROW(uWriter writer("test.txt", vs, "BED6"), no_fields_names);
	ASSERT_THROW(uWriter writer("test.txt", vs, "CUSTOM"), no_fields_names);
}

/*
 * Tests for uWriter ostream custom constructor: 
 *		uWriter(std::ostream* os, const std::vector<std::string>& fieldsNames, const std::string& type);
 *	Valid case:
 *		NormalUsage
 *	Invalid case:
 *		InvalidOstream	
 *		EmptyFieldsNames
 */

TEST(uWriterCustomOstreamConstructor, NormalUsage) {
	ostringstream oss;
	vector<string> vs;
	vs.push_back("CHR");
	ASSERT_NO_THROW(uWriter writer(&oss, vs, "BED4"));
	ASSERT_NO_THROW(uWriter writer(&oss, vs, "BED6"));
	ASSERT_NO_THROW(uWriter writer(&oss, vs, "CUSTOM"));
}

TEST(uWriterCustomOstreamConstructor, InvalidOstream) {
	ostringstream* oss = nullptr;
	vector<string> vs;
	vs.push_back("CHR");
	ASSERT_THROW(uWriter writer(oss, vs, "BED4"), std::runtime_error);
	ASSERT_THROW(uWriter writer(oss, vs, "BED6"), std::runtime_error);
	ASSERT_THROW(uWriter writer(oss, vs, "CUSTOM"), std::runtime_error);
}

TEST(uWriterCustomOstreamConstructor, EmptyFieldsNames) {
	ostringstream oss;
	vector<string> vs;
	ASSERT_THROW(uWriter writer(&oss, vs, "BED4"), no_fields_names);
	ASSERT_THROW(uWriter writer(&oss, vs, "BED6"), no_fields_names);
	ASSERT_THROW(uWriter writer(&oss, vs, "CUSTOM"), no_fields_names);
}

/*
 * Tests for the function:
 * 		void printString(const std::string& str);
 *	ValidCases:
 *		ValidString
 *		EmptyString
 */

TEST(uWriterPrintString, ValidString) {
	ostringstream oss;
	uWriter writer(&oss, "BED4");
	writer.printString("test string\n");
	ASSERT_EQ(oss.str(), "test string\n");
}

TEST(uWriterPrintString, EmptyString) {
	ostringstream oss;
	uWriter writer(&oss, "BED4");
	writer.printString("");
	ASSERT_EQ(oss.str(), "");
}
