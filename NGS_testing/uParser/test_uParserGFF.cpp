#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;

TEST(parserGFFTest, 7ITEMCOUNT)
{
    const int ITEMCOUNT=7;
    /**< Validate opening and first value*/
    try
    {
        uParser ourParser("../data/GFF/validsample.gff", "GFF");
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
        ASSERT_NO_THROW(uParser ourParser("../data/GFF/validsample.gff", "GFF"));
        uParser ourParser("../data/GFF/validsample.gff", "GFF");
        try {
            auto token=ourParser.getNextEntry();
        }
        catch(ugene_exception_base &e)
        {
            cout << fetchStringError(e);
        }
        catch(std::exception &e)
        {
            cout << e.what();
        }
    }
    uParser ourParser("../data/GFF/validsample.gff", "GFF");
    int count=0;
    while (ourParser.eof()==false)
    {
        uToken tempToken=ourParser.getNextEntry();
        count++;
    }
    EXPECT_EQ(count, ITEMCOUNT);
}


TEST(parserGFFTest, GTFSAMPLE)
{
    const int ITEMCOUNT=7;
    /**< Validate opening and first value*/
    try
    {
        uParser ourParser("../data/GTF/validsample.gtf", "GFF");
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
        ASSERT_NO_THROW(  uParser ourParser("../data/GTF/validsample.gtf", "GFF"));
        uParser ourParser("../data/GFF/validsample.gff", "GFF");
        try {
            auto token=ourParser.getNextEntry();
        }
        catch(ugene_exception_base &e)
        {
            cout << fetchStringError(e);
        }
        catch(std::exception &e)
        {
            cout << e.what();
        }
    }
    uParser ourParser("../data/GFF/validsample.gff", "GFF");
    int count=0;
    while (ourParser.eof()==false)
    {
        uToken tempToken=ourParser.getNextEntry();
        count++;
    }
    EXPECT_EQ(count, ITEMCOUNT);
}
