
#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;

TEST(uParserBedGraph_Ctr, openFile) {
    EXPECT_NO_THROW(uParser ourParser("../data/bedgraph/Small10K.bedgraph", "BEDGRAPH"));
}

TEST(uParserBedGraph_getNextEntry, fileCount) {
    const int numTokens=1000;
    uParser ourParser("../data/bedgraph/Small10K.bedgraph", "BEDGRAPH");
    int countTok=0;
    while (ourParser.eof()==false){
        ourParser.getNextEntry();
        countTok ++;
    }
    EXPECT_EQ(countTok,numTokens);
}

TEST(uParserBedGraph_getNextEntry, CommentENDL) {
    uParser ourParser("../data/bedgraph/smallCommentBrowser.bedgraph", "BEDGRAPH");
      int countTok=0;
    while (ourParser.eof()==false){
        EXPECT_NO_THROW(ourParser.getNextEntry());
        countTok ++;
    }
     EXPECT_EQ(countTok,6);
}


