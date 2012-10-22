
#include <gtest/gtest.h>
#include <string>
#include <sstream>

//TODO : Split tests by type of file??
/**< To test private member directly, changed private for public */

#include "IO/Parser/uParser.h"

using namespace std;
using namespace NGS;

TEST(parserHeaderSamTests, ValidHeader) {
    EXPECT_NO_THROW(uParser ourParser("../testFiles/SAM/valid_Header.sam", "SAM"));
}

TEST(parserHeaderSamTests, InvalidLN_VALUE) {
	EXPECT_THROW(uParser ourParser("../testFiles/SAM/invalid_LN.sam", "SAM"),invalid_header_param_throw); }

TEST(parserHeaderSamTests, MissingVN) {
	  EXPECT_THROW(uParser ourParser("../testFiles/SAM/missing_HD_VN.sam", "SAM"),uParser_invalid_Sam_header); }

TEST(parserHeaderSamTests, MissingSN) {
	EXPECT_THROW(uParser ourParser("../testFiles/SAM/missing_SQ_SN.sam", "SAM"),uParser_invalid_Sam_header); }

TEST(parserHeaderSamTests, InvalidLine) {
	 EXPECT_THROW(uParser ourParser("../testFiles/SAM/invalid_header.sam", "SAM"),uParser_invalid_Sam_header  ); }

TEST(parserHeaderSamTests, InvalidFilePath) {
	 EXPECT_THROW(uParser ourParser("../testFiles/SAM/nothinghere.sam", "SAM"),std::runtime_error  ); }

TEST(parserHeaderSamTests, MissingSQ) {
	EXPECT_THROW(uParser ourParser("../testFiles/SAM/missing_SQ_SN.sam", "SAM"),uParser_invalid_Sam_header); }

TEST(parserHeaderSamTests, getEntryUnitary) {

    const std::string CHR="chr21";
    const int START = 42653323;
    const int END=START+158-1;
    EXPECT_NO_THROW(uParser ourParser("../testFiles/SAM/UnitaryValid.sam", "SAM"));
	uParser ourParser("../testFiles/SAM/UnitaryValid.sam", "SAM");

    uToken Token =ourParser.getNextEntry();

    EXPECT_EQ(CHR,Token.getParam(token_param::CHR));
    EXPECT_EQ(START,std::stoi(Token.getParam(token_param::START_POS)));
    EXPECT_EQ(END,std::stoi(Token.getParam(token_param::END_POS)));


    //cout << Token.getParam(token_param::CHR) <<Token.getParam(token_param::START_POS) << Token.getParam(token_param::END_POS) <<endl;
}

	TEST(parserHeaderSamTests, 5ITEMCOUNT) {
    const int ITEMCOUNT=5;
    /**< Validate opening and first value*/
{
    EXPECT_NO_THROW(uParser ourParser("../testFiles/SAM/fiveCountValid.sam", "SAM"));
    uParser ourParser("../testFiles/SAM/fiveCountValid.sam", "SAM");
    EXPECT_NO_THROW(auto token=ourParser.getNextEntry());
}
	uParser ourParser("../testFiles/SAM/fiveCountValid.sam", "SAM");

    int count=0;
    while (ourParser.eof()==false){
        uToken tempToken=ourParser.getNextEntry();
        count++;
    }
    EXPECT_EQ(count, ITEMCOUNT);
    }
