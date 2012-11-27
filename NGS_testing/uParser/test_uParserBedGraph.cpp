
#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;

TEST(parserOpenBedgraph, openFile) {
    EXPECT_NO_THROW(uParser ourParser("../data/bedgraph/Small10K.bedgraph", "BEDGRAPH"));
}


TEST(parserOpenBedgraph, fileCount) {
    const int numTokens=1000;
    uParser ourParser("../data/bedgraph/Small10K.bedgraph", "BEDGRAPH");
    int countTok=0;
    while (ourParser.eof()==false){
        ourParser.getNextEntry();
        countTok ++;
    }
    EXPECT_EQ(countTok,numTokens);

}
