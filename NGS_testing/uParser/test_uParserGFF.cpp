#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;



TEST(parserGFF_CTR, VALIDPATH)
{
    EXPECT_NO_THROW(uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF"));

}

TEST(parserGFF_CTR, INVALIDPATH)
{
    EXPECT_ANY_THROW(uParser ourParser("../data/GFF/nothinghere.gff", "UCSCGFF"));
}

TEST(parserGFF_getNextEntry, ENDLD)
{
    uParser ourParser("../data/GFF/endl.gff", "UCSCGFF");
    while (ourParser.eof()==false)
    {
       ASSERT_NO_THROW( uToken tempToken=ourParser.getNextEntry());
    }
}

TEST(parserGFF_getNextEntry, UCSCHEADER)
{
    uParser ourParser("../data/GFF/UCSCSample.gff", "UCSCGFF");
    while (ourParser.eof()==false)
    {
       ASSERT_NO_THROW( uToken tempToken=ourParser.getNextEntry());
    }
}


TEST(parserGFF, 7ITEMCOUNT)
{
    const int ITEMCOUNT=7;
    /**< Validate opening and first value*/
    try
    {
        uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF");
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
        ASSERT_NO_THROW(uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF"));
        uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF");
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
    uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF");
    int count=0;
    while (ourParser.eof()==false)
    {
        uToken tempToken=ourParser.getNextEntry();
        count++;
    }
    EXPECT_EQ(count, ITEMCOUNT);
}

TEST(parserGFF, GFFSAMPLE)
{
    const int ITEMCOUNT=7;
    /**< Validate opening and first value*/
    try
    {
        uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF");
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
        ASSERT_NO_THROW(  uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF"));
        uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF");
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
    uParser ourParser("../data/GFF/validsample.gff", "UCSCGFF");
    int count=0;
    while (ourParser.eof()==false)
    {
        uToken tempToken=ourParser.getNextEntry();
        count++;
    }
    EXPECT_EQ(ITEMCOUNT,count);
}
