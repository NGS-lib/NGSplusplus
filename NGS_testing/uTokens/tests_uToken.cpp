#include <gtest/gtest.h>
#include "NGS++.h"
#include <string>
#include <sstream>

using namespace std;
using namespace NGS;

// TODO: test CIGAR containing only a '*'
/*
 * Constructor tests:
 *		uToken(std::istream& paramList);
 * 	Valid cases:
 *		MandatoryArgumentsOnly
 *		STRAND_IsValid
 *		MAP_SCORE_IsValid_lowerBound
 *		MAP_SCORE_IsValid_higherBound
 *		START_POS_equal_END_POS
 *		PHRED_SCORE_SEQUENCE_AreValid
 *		CIGAR_SEQUENCE_AreValid
 *		FLAGS_IsValid_lowerBound
 *		FLAGS_IsValid_higherBound
 *		AllArgumentsAreValids
 *	Invalid cases:
 * 		EmptyArgument
 *		Invalid_param_token
 *		Invalid_param_token
 *		MissingMandatoryArgument_CHR
 *		MissingMandatoryArgument_START_POS
 *		MissingMandatoryArgument_END_POS
 *		Invalid_START_POS
 *		Invalid_END_POS
 *		Invalid_STRAND
 *		Invalid_MAP_SCORE_lowerBound
 *		Invalid_MAP_SCORE_higherBound
 *		Invalid_FLAGS_lowerBound
 *		Invalid_FLAGS_higherBound
 *		Invalid_SEQUENCE
 *		Invalid_CIGAR_2alphabeticCharInARow
 *		Invalid_CIGAR_EndWithNumber
 *		Invalid_START_END_Combination
 *		SEQUENCE_PHRED_lengthDontMatch
 *		SEQUENCE_CIGAR_lengthDontMatch
 */

TEST(uTokenConstructor, MandatoryArgumentsOnly) {
	stringstream ss;
	ss << "START_POS\t1\n" << "END_POS\t21\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, STRAND_IsValid) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, MAP_SCORE_IsValid_lowerBound) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "MAP_SCORE\t0\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, MAP_SCORE_IsValid_higherBound) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "MAP_SCORE\t255\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, START_POS_equal_END_POS) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t1\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, PHRED_SCORE_SEQUENCE_AreValid) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "PHRED_SCORE\t###########\n" << "SEQUENCE\tACGTN.acgtn\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, CIGAR_SEQUENCE_AreValid) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, FLAGS_IsValid_lowerBound) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "FLAGS\t1\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, FLAGS_IsValid_higherBound) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "FLAGS\t65535\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, AllArgumentsAreValids) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, EmptyArgument) {
	stringstream ss;
	ss << "CHR\t\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_param_token) {
	stringstream ss;
	ss << "CH\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ASSERT_THROW(uToken Token(ss), uToken_exception_base);
}

TEST(uTokenConstructor, Missing_CHR) {
	stringstream ss;
	ss << "START_POS\t1\n" << "END_POS\t21\n";
	ASSERT_NO_THROW(uToken Token(ss));
}

TEST(uTokenConstructor, MissingMandatoryArgument_START_POS) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "END_POS\t21\n";
	ASSERT_THROW(uToken Token(ss), invalid_uToken_throw);
}

TEST(uTokenConstructor, MissingMandatoryArgument_END_POS) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n";
	ASSERT_THROW(uToken Token(ss), invalid_uToken_throw);
}

TEST(uTokenConstructor, Invalid_START_POS) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t-1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_END_POS) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t-1\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_STRAND) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t=\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_MAP_SCORE_lowerBound) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t-1\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_MAP_SCORE_higherBound) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t256\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_FLAGS_lowerBound) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t-1\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_FLAGS_higherBound) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t65536\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_SEQUENCE) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCx\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_CIGAR_2alphabeticCharInARow) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=1XX\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_CIGAR_EndWithNumber) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=123\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_value_throw);
}

TEST(uTokenConstructor, Invalid_START_END_Combination) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t21\n" << "END_POS\t1\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_uToken_throw);
}

TEST(uTokenConstructor, SEQUENCE_PHRED_lengthDontMatch) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t###################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_uToken_throw);
}

TEST(uTokenConstructor, SEQUENCE_CIGAR_lengthDontMatch) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t1M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	ASSERT_THROW(uToken Token(ss), invalid_uToken_throw);
}

/*
 * Test for the function:
 *		std::string getParam(token_param name);
 *	Valid case:
 *		ArgumentExists
 *	Invalid case:
 *		ArgumentDoesNotExist
 */

 TEST(uTokenGetParam, ArgumentExists) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	uToken Token(ss);
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "21");
	ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
	ASSERT_EQ(Token.getParam(token_param::MAP_SCORE), "255");
	ASSERT_EQ(Token.getParam(token_param::PHRED_SCORE), "####################");
	ASSERT_EQ(Token.getParam(token_param::CIGAR), "2M3I2X1=12X");
	ASSERT_EQ(Token.getParam(token_param::SEQUENCE), "ACGTN.acgtn.ACGTGTCN");
	ASSERT_EQ(Token.getParam(token_param::SEQ_NAME), "ab00001");
	ASSERT_EQ(Token.getParam(token_param::FLAGS), "256");
}

TEST(uTokenGetParam, ArgumentDoesNotExist) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	uToken Token(ss);
	ASSERT_EQ(Token.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token.getParam(token_param::START_POS), "1");
	ASSERT_EQ(Token.getParam(token_param::END_POS), "21");
	ASSERT_THROW(Token.getParam(token_param::STRAND),param_not_found);
	ASSERT_THROW(Token.getParam(token_param::MAP_SCORE), param_not_found);
	ASSERT_THROW(Token.getParam(token_param::PHRED_SCORE), param_not_found);
	ASSERT_THROW(Token.getParam(token_param::CIGAR), param_not_found);
	ASSERT_THROW(Token.getParam(token_param::SEQUENCE),param_not_found);
	ASSERT_THROW(Token.getParam(token_param::SEQ_NAME), param_not_found);
	ASSERT_THROW(Token.getParam(token_param::FLAGS), param_not_found);
}

/*
 * Test for the function:
 *		std::string getParam(const std::string& name) const;
 *	Valid cases:
 *		ArgumentExistsTokenParam
 *		ArgumentExistsCustomParam
 *	Invalid case:
 *		ArgumentDoesNotExist
 */


TEST(uTokenGetParamString, ArgumentExistsTokenParam) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	uToken Token(ss);
	ASSERT_EQ(Token.getParam("CHR"), "chr1");
	ASSERT_EQ(Token.getParam("START_POS"), "1");
	ASSERT_EQ(Token.getParam("END_POS"), "21");
	ASSERT_EQ(Token.getParam("STRAND"), "+");
	ASSERT_EQ(Token.getParam("MAP_SCORE"), "255");
	ASSERT_EQ(Token.getParam("PHRED_SCORE"), "####################");
	ASSERT_EQ(Token.getParam("CIGAR"), "2M3I2X1=12X");
	ASSERT_EQ(Token.getParam("SEQUENCE"), "ACGTN.acgtn.ACGTGTCN");
	ASSERT_EQ(Token.getParam("SEQ_NAME"), "ab00001");
	ASSERT_EQ(Token.getParam("FLAGS"), "256");
}

TEST(uTokenGetParamString, ArgumentExistsCustomParam) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "JUNK\tabcd\n";
	uToken Token(ss,true);
	ASSERT_EQ(Token.getParam("CHR"), "chr1");
	ASSERT_EQ(Token.getParam("START_POS"), "1");
	ASSERT_EQ(Token.getParam("END_POS"), "21");
	ASSERT_EQ(Token.getParam("JUNK"), "abcd");
}

TEST(uTokenGetParamString, ArgumentDoesNotExist) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	uToken Token(ss);
	ASSERT_EQ(Token.getParam("CHR"), "chr1");
	ASSERT_EQ(Token.getParam("START_POS"), "1");
	ASSERT_EQ(Token.getParam("END_POS"), "21");
	ASSERT_THROW(Token.getParam("STRAND"),param_not_found);
	ASSERT_THROW(Token.getParam("MAP_SCORE"), param_not_found);
	ASSERT_THROW(Token.getParam("PHRED_SCORE"), param_not_found);
	ASSERT_THROW(Token.getParam("CIGAR"), param_not_found);
	ASSERT_THROW(Token.getParam("SEQUENCE"),param_not_found);
	ASSERT_THROW(Token.getParam("SEQ_NAME"), param_not_found);
	ASSERT_THROW(Token.getParam("FLAGS"), param_not_found);
}

TEST(uToken_paramCount, ONE) {

    stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	uToken Token(ss);
	EXPECT_EQ(1, Token.paramCount(token_param::START_POS));
}

TEST(uToken_paramCount, Several) {

    stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss <<"START_POS\t1\n" << "END_POS\t21\n";
	ss <<"START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	uToken Token(ss);
	EXPECT_EQ(3, Token.paramCount(token_param::START_POS));
}
TEST(uToken_paramCount, None) {

    stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	uToken Token(ss);
	EXPECT_EQ(0, Token.paramCount(token_param::CUST_3));
}
/*
 * Test for the checking if param is Set:
 *		bool uToken::_isParamSet;
 */

TEST(uTokenGetParam, ArgumentIsSet) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	uToken Token(ss);
	ASSERT_TRUE(Token.isParamSet(token_param::CHR));
	ASSERT_TRUE(Token.isParamSet(token_param::START_POS));
	ASSERT_TRUE(Token.isParamSet(token_param::END_POS));
}
TEST(uTokenGetParam, ArgumentisNotSet) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	uToken Token(ss);
	ASSERT_FALSE(Token.isParamSet(token_param::CIGAR));
	ASSERT_FALSE(Token.isParamSet(token_param::SEQUENCE));
}

/*
 * Test for the equality operator overloading:
 *		uToken& operator=(uToken const& assign_from);
 *	Valid cases:
 *		UsedAfterDeclaration
 *		UsedDuringDeclaration
 */

TEST(uTokenEqualityOperatorOverloading, UsedAfterDeclaration) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	uToken Token(ss);
	ss.clear();
	ss << "CHR\tchr2\n" << "START_POS\t11\n" << "END_POS\t31\n";
	uToken Token2(ss);
	ASSERT_EQ(Token2.getParam(token_param::CHR), "chr2");
	ASSERT_EQ(Token2.getParam(token_param::START_POS), "11");
	ASSERT_EQ(Token2.getParam(token_param::END_POS), "31");
	Token2 = Token;
	ASSERT_EQ(Token2.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token2.getParam(token_param::START_POS), "1");
	ASSERT_EQ(Token2.getParam(token_param::END_POS), "21");
	ASSERT_EQ(Token2.getParam(token_param::STRAND), "+");
	ASSERT_EQ(Token2.getParam(token_param::MAP_SCORE), "255");
	ASSERT_EQ(Token2.getParam(token_param::PHRED_SCORE), "####################");
	ASSERT_EQ(Token2.getParam(token_param::CIGAR), "2M3I2X1=12X");
	ASSERT_EQ(Token2.getParam(token_param::SEQUENCE), "ACGTN.acgtn.ACGTGTCN");
	ASSERT_EQ(Token2.getParam(token_param::SEQ_NAME), "ab00001");
	ASSERT_EQ(Token2.getParam(token_param::FLAGS), "256");
}

TEST(uTokenEqualityOperatorOverloading, UsedDuringDeclaration) {
	stringstream ss;
	ss << "CHR\tchr1\n" << "START_POS\t1\n" << "END_POS\t21\n";
	ss << "STRAND\t+\n" << "MAP_SCORE\t255\n" << "PHRED_SCORE\t####################\n";
	ss << "CIGAR\t2M3I2X1=12X\n" << "SEQUENCE\tACGTN.acgtn.ACGTGTCN\n";
	ss << "SEQ_NAME\tab00001\n" << "FLAGS\t256\n";
	uToken Token(ss);
	uToken Token2 = Token;
	ASSERT_EQ(Token2.getParam(token_param::CHR), "chr1");
	ASSERT_EQ(Token2.getParam(token_param::START_POS), "1");
	ASSERT_EQ(Token2.getParam(token_param::END_POS), "21");
	ASSERT_EQ(Token2.getParam(token_param::STRAND), "+");
	ASSERT_EQ(Token2.getParam(token_param::MAP_SCORE), "255");
	ASSERT_EQ(Token2.getParam(token_param::PHRED_SCORE), "####################");
	ASSERT_EQ(Token2.getParam(token_param::CIGAR), "2M3I2X1=12X");
	ASSERT_EQ(Token2.getParam(token_param::SEQUENCE), "ACGTN.acgtn.ACGTGTCN");
	ASSERT_EQ(Token2.getParam(token_param::SEQ_NAME), "ab00001");
	ASSERT_EQ(Token2.getParam(token_param::FLAGS), "256");
}
