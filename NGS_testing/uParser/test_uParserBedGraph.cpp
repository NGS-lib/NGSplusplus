
#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;

const string SMALL1K="../data/bedgraph/Small10K.bedgraph";
const string UCSCBROWSER="../data/bedgraph/smallCommentBrowser.bedgraph";

TEST(uParserBedGraph_Ctr, openFile) {
    EXPECT_NO_THROW(uParser ourParser(SMALL1K, "BEDGRAPH"));
}

TEST(uParserBedGraph_getNextEntry, fileCount) {
    const int numTokens=1000;
    uParser ourParser(SMALL1K, "BEDGRAPH");
    int countTok=0;
    while (ourParser.eof()==false){
        ourParser.getNextEntry();
        countTok++;
    }
    EXPECT_EQ(countTok,numTokens);
}

TEST(uParserBedGraph_getNextEntry, CommentENDL) {
    uParser ourParser(UCSCBROWSER, "BEDGRAPH");
      int countTok=0;
    while (ourParser.eof()==false){
        EXPECT_NO_THROW(ourParser.getNextEntry());
        countTok ++;
    }
     EXPECT_EQ(countTok,6);
}

TEST(uParserBedGraph_getNextEntry, MAKEBASIC) {
    uParser ourParser(SMALL1K, "BEDGRAPH");
    int countTok=0;

    vector<long long> startVec={10301,10321,10326,10331,10336,10345,10348,10357,10362,10367};
    vector<long long> endVec={10321,10326,10331,10336,10345,10348,10357,10362,10367,10371};
    vector<float> scoreVec={4.0,4.2,4.2,4.3,4.4,4.4,4.5,4.3,4.2,4.1};
    /**< Check first ten */
    vector<uBasicNGS> tenBasic;
    while (ourParser.eof()==false){
        tenBasic.push_back(uBasicNGS(ourParser.getNextEntry()));
        countTok ++;
        if (countTok==10)
            break;
    }

    for(int i=0;i<10;i++)
    {
        ASSERT_EQ(tenBasic.at(0).getStart(),startVec.at(0));
        ASSERT_EQ(tenBasic.at(0).getEnd(),endVec.at(0));
        ASSERT_EQ(tenBasic.at(0).getScore(),scoreVec.at(0));
    }

}
