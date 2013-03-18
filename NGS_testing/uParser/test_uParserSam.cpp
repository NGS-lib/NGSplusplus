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

TEST(uParserSam_getNextEntry, MAKEBASICTEN) {
    uParser ourParser(BIG_FILE, "SAM");
    int countTok=0;

    vector<long long> startVec={56262,61780,64367,67342,71763,71788,71788,71789,71792,72128};
    vector<long long> endVec={56262+39,61780+39,64367+39,67342+39,71763+39,71788+39,71788+39,71789+39,71792+39,72128+39};
    vector<string> seqVec={
    "CAAGGACCGGACGGATTCACAGCAGAATTCTACCAGACAT",
    "TTCTTTTTAAACCAGCCGGCTGTGGCGGCATATGCCTGTA",
    "AGCAGTCCCTGATCCTGTTTCTCCTGACTGGACAAGACCT",
    "AGTCTCAGATCTTTGGGAGGCAAAGGTGGGCAGATCACCT",
    "TCGTACGCCCATGGAGTCTCGCTGATTGCTAGCACAACCG",
    "TTGCTAGCACAACCGTCTGAGATCAAACTGCAAGGCGGCA",
    "TTGCTAGCACAACCGTCTGAGATCAAACTGCAAGGCGGCA",
    "TGCTAGCACAACCGTCTGAGATCAAACTGCAAGGCGGCAG",
    "TAGCACAACCGTCTGAGATCAAACTGCAAGGCGGCAGCGA",
    "TGGGAGGCACCCCCCAGCAGGGCACACTGACACCTCACAT",
    "ACCCCCCAGCAGGGCACACTGACACCTCACATGGCAGGGT"};

    vector<string> seqNameVec={
    "HWI-ST333_0095_FC:5:1102:6572:40808#GTACGT",
    "HWI-ST333_0095_FC:3:1207:10472:200099#GGTACC",
    "HWI-ST333_0095_FC:3:1208:11369:155611#GGTACC",
    "HWI-ST333_0095_FC:5:2104:14545:176320#GTACGT",
    "HWI-ST333_0095_FC:5:1107:10883:5991#GTACGT",
    "HWI-ST333_0095_FC:3:1103:18911:137521#GGTACC",
    "HWI-ST333_0095_FC:5:1106:9536:83045#GTACGT",
    "HWI-ST333_0095_FC:5:1207:3471:37381#GTACGT",
    "HWI-ST333_0095_FC:3:2205:17931:44824#GGTACC",
    "HWI-ST333_0095_FC:5:1104:1435:126724#GTACGT"};

    vector<string> phred={
    "BBBBBBBBBBB^UYUcc`^`bb_b\\[RUXPccaccccccc",
    "gdgeddddffcggefggfegggdgggggggggggg^ffff",
    "gggggggggggggggggggggggggggggggggggggggg",
    "gggggggggdggggggggggfggggd\\eecb`cbbgdfgg",
    "ggggeggggegggedgggegdggggfgggggggggfefef",
    "gggggggggggggggggggggeggggggggggegggggge",
    "eggggegfdggdgggggggggggggggggggggggggggg",
    "ggggggggggggfggggggggefgfggggggggggggege",
    "egggeggegggegggggggggggggggggggggggggggg",
    "gggggggfggggggggggggggggfggeggggggggggge"
    };

    /**< Check first ten */
    vector<uTags> tenTags;
    while (ourParser.eof()==false){
        tenTags.push_back(uTags(ourParser.getNextEntry()));
        countTok ++;
        if (countTok==10)
            break;
    }

    for(int i=0;i<10;i++)
    {
        ASSERT_EQ(startVec.at(i),tenTags.at(i).getStart());
        ASSERT_EQ(endVec.at(i),tenTags.at(i).getEnd());
        ASSERT_EQ(seqVec.at(i),tenTags.at(i).getSequence());
        ASSERT_EQ(seqNameVec.at(i),tenTags.at(i).getName());
        ASSERT_EQ(phred.at(i),tenTags.at(i).getPhred());

    }

}


TEST(parserSamTests, BIGPARSE30K)
{
    EXPECT_NO_THROW(uParser ourParser(BIG_FILE, "SAM");
    int count=30000;
    int cur=0;
    while (cur!=count)
        {
            uToken tempToken=ourParser.getNextEntry();
            cur++;
        }
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
