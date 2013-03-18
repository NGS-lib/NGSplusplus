#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;
using namespace BamTools;
const string BAMREAL="../data/BAM/H2AZ.bam";
vector<int> startPosBam={56262,61780,64367,67342,71763};
TEST(uParserBam_getNextEntry, REALDATA)
{
   uParser ourParser(BAMREAL, "BAM");
   for(int i=0; i<10; i++)
   {
      ASSERT_NO_THROW(ourParser.getNextEntry());
   }
}

TEST(uParserBam_getNextEntry, COMPARETOKENTAG)
{
   uParser ourParser(BAMREAL, "BAM");

   for(int i=0; i<5; i++)
   {
       uToken curToken = ourParser.getNextEntry();
       ASSERT_EQ(startPosBam.at(i), uTags(curToken).getStart());
       ASSERT_EQ(uTags(curToken).getStart(),std::stoi(curToken.getParam(token_param::START_POS)));
   }
}

