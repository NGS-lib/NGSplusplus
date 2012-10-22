
#include <gtest/gtest.h>
#include <string>
#include <sstream>

//TODO : Split tests by type of file??
/**< To test private member directly, changed private for public */
//#define private public
//#include "IO/Parser/uParser.h" // TODO: When the new parser will be functional, we will change it
#include "IO/Parser/uParser.h"

using namespace std;
using namespace NGS;

TEST(newParserWig, CorrectlyFormatedVariableWIG) {

	uParser ourParser("../wig/correctVariable.wig", "WIG");
	int count=0;
	uToken Token = ourParser.getNextEntry();
	count++;
	EXPECT_EQ(Token.getParam(token_param::CHR), "chr19");
	EXPECT_EQ(Token.getParam(token_param::START_POS), "49304701");
	EXPECT_EQ(Token.getParam(token_param::END_POS), std::to_string(49304701+150));
	EXPECT_EQ(Token.getParam(token_param::SCORE), "10");
	Token = ourParser.getNextEntry();
	count++;
	EXPECT_EQ(Token.getParam(token_param::CHR), "chr19");
	EXPECT_EQ(Token.getParam(token_param::START_POS), "49304901");
	EXPECT_EQ(Token.getParam(token_param::END_POS), std::to_string(49304901+150));
	EXPECT_EQ(Token.getParam(token_param::SCORE), "12.5");
	//ASSERT_EQ(Token.getParam(token_param::STRAND), "+");
EXPECT_NO_THROW(
	while(!(ourParser.eof())){
    Token = ourParser.getNextEntry();
	count++;
	}
);
EXPECT_EQ(count, 9);
}
TEST(newParserWig, CorrectlyFormatedFixedWIG) {
    int count=0;
    uParser ourParser("../wig/correctFixed.wig", "WIG");
    const int START=49307401;
	const int STEP=300;
	const int SPAN=200;
     uToken Token = ourParser.getNextEntry();

	count++;
	EXPECT_EQ(Token.getParam(token_param::CHR), "chr19");
	EXPECT_EQ( std::to_string(START),Token.getParam(token_param::START_POS));
	EXPECT_EQ(Token.getParam(token_param::END_POS), std::to_string(START+SPAN));
	EXPECT_EQ(Token.getParam(token_param::SCORE), std::to_string(1000));
	Token = ourParser.getNextEntry();
	count++;
	EXPECT_EQ(Token.getParam(token_param::CHR), "chr19");
	EXPECT_EQ(Token.getParam(token_param::START_POS), std::to_string(START+STEP));
	EXPECT_EQ(Token.getParam(token_param::END_POS), std::to_string(START+STEP+SPAN));
	EXPECT_EQ(Token.getParam(token_param::SCORE), std::to_string(900));
    count++;
    Token = ourParser.getNextEntry();
	EXPECT_EQ(Token.getParam(token_param::CHR), "chr19");
	EXPECT_EQ(Token.getParam(token_param::START_POS), std::to_string(START+STEP*2));
	EXPECT_EQ(Token.getParam(token_param::END_POS), std::to_string( START+STEP*2+SPAN));
	EXPECT_EQ(Token.getParam(token_param::SCORE), std::to_string(800));
EXPECT_NO_THROW(
	while(!(ourParser.eof())){
    Token = ourParser.getNextEntry();

	count++;
	}
);
EXPECT_EQ(count, 10);
}

TEST(newParserWig, CorrectFormatBoth) {
    int count=0;
    uParser ourParser("../wig/correct.wig", "WIG");
    const int STEP=300;
    vector<int> compareStart{49304701,49304901,49305401,49305601,49305901,49306081,49306301,49306691,49307871,49307401,49307401+STEP,49307401+STEP*2,49307401+STEP*3,49307401+STEP*4,49307401+STEP*5,49307401+STEP*6,49307401+STEP*7,49307401+STEP*8,49307401+STEP*9};
   // vector<int> compareEnd{}
    vector<float> comparescores{10.0,12.5,15.0,17.5,20.0,17.5,15.0,12.5,10.0,1000,900,800,700,600,500,400,300,200,100};
    vector<float> scores;
    vector<int> starts;
EXPECT_NO_THROW(
	while(!(ourParser.eof())){
    uToken Token = ourParser.getNextEntry();
	scores.push_back(std::stof(Token.getParam(token_param::SCORE)));
	starts.push_back(std::stoi(Token.getParam(token_param::START_POS)));
	count++;
	}

);
EXPECT_EQ(count, 19);
EXPECT_EQ(comparescores, scores);
EXPECT_EQ(compareStart, starts);
}
TEST(newParserWig, noChromWig) {

EXPECT_ANY_THROW(
     uParser ourParser("../wig/nochromFixed.wig", "WIG");
);

}
TEST(newParserWig, noStartFixed) {

EXPECT_ANY_THROW(
    uParser ourParser("../wig/noStartFixed.wig", "WIG");
);

}

TEST(newParserWig, noDefinition) {

EXPECT_ANY_THROW(
    uParser ourParser("../wig/no_definition.wig", "WIG");
);

}

TEST(newParserWig, incorrectLines) {
    uParser ourParser("../wig/incorrect.wig", "WIG");
    EXPECT_ANY_THROW(

    while(!(ourParser.eof())){
         uToken Token = ourParser.getNextEntry();
    }
);
}