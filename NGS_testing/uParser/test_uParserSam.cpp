#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;

const string VALIDHEADER="../data/SAM/valid_Header.sam";
const string INVALIDLN="../data/SAM/invalid_LN.sam";
const string MISSINGHD="../data/SAM/missing_HD_VN.sam";
const string MISSINGSQ="../data/SAM/missing_SQ_SN.sam";
const string INVALID_HEAD="../data/SAM/invalid_header.sam";
const string BIG_FILE="../data/SAM/H2AZFinal_downsample.sam";
const string UNITARY_VALID="../data/SAM/UnitaryValid.sam";
const string FIVE_COUNT="../data/SAM/fiveCountValid.sam";

TEST(parserHeaderSamTests, ValidHeader)
{
    EXPECT_NO_THROW(uParser ourParser(VALIDHEADER, "SAM"));
}

TEST(parserHeaderSamTests, InvalidLN_VALUE)
{
    EXPECT_THROW(uParser ourParser(INVALIDLN, "SAM"),invalid_header_param_throw);
}

TEST(parserHeaderSamTests, MissingVN)
{
    EXPECT_THROW(uParser ourParser(MISSINGHD, "SAM"),uParser_invalid_Sam_header);
}

TEST(parserHeaderSamTests, MissingSN)
{
    EXPECT_THROW(uParser ourParser(MISSINGSQ, "SAM"),uParser_invalid_Sam_header);
}

TEST(parserHeaderSamTests, InvalidLine)
{
    EXPECT_THROW(uParser ourParser(INVALID_HEAD, "SAM"),uParser_invalid_Sam_header  );
}

TEST(parserHeaderSamTests, InvalidFilePath)
{
    EXPECT_THROW(uParser ourParser("../data/SAM/nothinghere.sam", "SAM"),std::runtime_error  );
}

TEST(parserHeaderSamTests, MissingSQ)
{
    EXPECT_THROW(uParser ourParser(MISSINGSQ, "SAM"),uParser_invalid_Sam_header);
}

TEST(parserHeaderSamTests, getEntryUnitary)
{

    const std::string CHR="chr21";
    const int START = 42653323;
    const int END=START+158-1;
    EXPECT_NO_THROW(uParser ourParser(UNITARY_VALID, "SAM"));
    uParser ourParser(UNITARY_VALID, "SAM");
    uToken Token =ourParser.getNextEntry();
    EXPECT_EQ(CHR,Token.getParam(token_param::CHR));
    EXPECT_EQ(START,std::stoi(Token.getParam(token_param::START_POS)));
    EXPECT_EQ(END,std::stoi(Token.getParam(token_param::END_POS)));


}


TEST(parserSamTests, BIGFILE)
{
    EXPECT_NO_THROW(uParser ourParser(BIG_FILE, "SAM");

    while (ourParser.eof()==false)
        {uToken tempToken=ourParser.getNextEntry(); }
    );

}


TEST(parserHeaderSamTests, 5ITEMCOUNT)
{
    const int ITEMCOUNT=5;
    /**< Validate opening and first value*/

    try
    {
        uParser ourParser(FIVE_COUNT, "SAM");
    }
    catch(ugene_exception_base &e)
    {
        cout << fetchStringError(e);
    }
    catch(std::exception &e)
    {
        cout << e.what();
    }

    {
        ASSERT_NO_THROW(uParser ourParser(FIVE_COUNT, "SAM"));
        uParser ourParser(FIVE_COUNT, "SAM");
        EXPECT_NO_THROW(auto token=ourParser.getNextEntry());
    }
    uParser ourParser(FIVE_COUNT, "SAM");

    int count=0;
    while (ourParser.eof()==false)
    {
        uToken tempToken=ourParser.getNextEntry();
        count++;
    }
    EXPECT_EQ(count, ITEMCOUNT);
}
