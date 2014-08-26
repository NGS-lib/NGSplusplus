#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;
const string BWReal="/home/local/USHERBROOKE/nora2001/NGSplusplus/NGS_testing/data/BW/mginput.bigwig";

TEST(uParserBW_getNextEntry, REALDATA)
{
   uParser ourParser(BWReal, "BW");

   while(ourParser.eof()==false)
   {
      ASSERT_NO_THROW(ourParser.getNextEntry());
   }
}

