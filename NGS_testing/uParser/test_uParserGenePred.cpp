#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "NGS++.h"

using namespace std;
using namespace NGS;

const string HG18UCSC="../data/GENEPRED/hg18GenePred.genePred";
vector<int> startPosHg18={1115,1115,4268,5658,6720};
vector<int> EndPosHg18={4121,4272,6628,7231,7614};

vector<int> exonStartFirst={1115,2475,3083};
vector<int> exonEndFirst={2090,2584,4121};

vector<int> exonStartSecond={1115,2475};
vector<int> exonEndSecond={2090,4272};



const string HG18UCSCENDL="../data/GENEPRED/hg18endl.genePred";

const string HG18NOHEADER="../data/GENEPRED/hg18endlnoheader.genePred";



TEST(uParserGenePred_ConstructorFilename, ValidFileName)
{
    ASSERT_NO_THROW(uParser Parser(HG18UCSC, "GENEPRED"));
}

TEST(uParserGenePred_getNextEntry, ParseOne)
{
    uParser Parser(HG18UCSC, "GENEPRED");
    auto token = Parser.getNextEntry();
}

TEST(uParserGenePred_getNextEntry, ENDL)
{
    uParser Parser(HG18UCSCENDL, "GENEPRED");
    uToken token;
    while (Parser.eof()==false){
      ASSERT_NO_THROW(token = Parser.getNextEntry(););
    }
}
TEST(uParserGenePred_getNextEntry, NOHEADER)
{
    uParser Parser(HG18NOHEADER, "GENEPRED");
        uToken token;
    while (Parser.eof()==false){
      ASSERT_NO_THROW(token = Parser.getNextEntry(););
    }
}



TEST(uParserGenePred_getNextEntry, Parse10K)
{
    uParser Parser(HG18UCSC, "GENEPRED");
    int loadCount=10000;
    int cur=0;
    while (cur!=loadCount)
    {
    ASSERT_NO_THROW(auto token = Parser.getNextEntry());
    cur++;
    }

}

TEST(uParserGenePred_getNextEntry, VALIDATE5)
{
    uParser Parser(HG18UCSC, "GENEPRED");
    for (int i=0; i<5; i++)
    {
         uToken token = Parser.getNextEntry();
         ASSERT_EQ(startPosHg18.at(i),std::stoi(token.getParam(token_param::START_POS)));
         ASSERT_EQ(EndPosHg18.at(i),std::stoi(token.getParam(token_param::END_POS)));
    }
}

TEST(uParserGenePred_getNextEntry, VALIDATEEXON2)
{
     uParser Parser(HG18UCSC, "GENEPRED");
     uToken token1 = Parser.getNextEntry();
     /**< Coding is first feature, we ignore */
     for (int i=2 ; (i<token1.paramCount(token_param::START_POS));i++ ){
        ASSERT_EQ(exonStartFirst.at(i-2),std::stoi(token1.getParam(token_param::START_POS,i)));
        ASSERT_EQ(exonEndFirst.at(i-2),std::stoi(token1.getParam(token_param::END_POS,i)));
     }

    uToken token2 = Parser.getNextEntry();
    for (int i=2 ; (i<token2.paramCount(token_param::START_POS));i++ ){
        ASSERT_EQ(exonStartSecond.at(i-2),std::stoi(token2.getParam(token_param::START_POS,i)));
        ASSERT_EQ(exonEndSecond.at(i-2),std::stoi(token2.getParam(token_param::END_POS,i)));
     }
}
