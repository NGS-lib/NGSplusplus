
#include <gtest/gtest.h>
#include <string>
#include <sstream>

//TODO : Split tests by type of file??
/**< To test private member directly, changed private for public */
//#define private public
//#include "IO/Parser/uParser.h" // TODO: When the new parser will be functional, we will change it
#include "IO/Parser/iParser.h"

using namespace std;
using namespace NGS;

TEST(newParserGetNext, CorrectlyFormatedVariableWIG) {
     std::cerr <<"Hello" <<std::endl;
	Parser ourParser("./wig/correctVariable.wig", "WIG");
    std::cerr <<"Hello before Parser" <<std::endl;
	int count=0;
	uToken Token = ourParser.getNextEntry();
	count++;
	EXPECT_EQ(Token.getParam(token_param::CHR), "chr19");
	EXPECT_EQ(Token.getParam(token_param::START_POS), "49304701");
	EXPECT_EQ(Token.getParam(token_param::END_POS), std::to_string(49304701+150));
	EXPECT_EQ(Token.getParam(token_param::SCORE), "10");
   std::cerr <<"Hello" <<std::endl;
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
std::cerr <<"Hello" <<std::endl;
}
