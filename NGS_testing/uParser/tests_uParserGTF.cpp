#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include "NGS++.h"

using namespace std;
using namespace NGS;

//Validator
//http://mblab.wustl.edu/software.html

const string UCSCKNOWNGENE="../data/GTF/knownGene.gtf";
const string ENDL="../data/GTF/endl.gtf";
const string HEADER="../data/GTF/UCSCheader.gtf";


TEST(parserGTF_CTR, VALIDPATH)
{
    /**< Validate opening and first value*/
    uParser ourParser(UCSCKNOWNGENE, "GTF");
}


TEST(parserGTF_CTR, INVALIDPATH)
{
    /**< Validate opening and first value*/
    EXPECT_ANY_THROW(uParser ourParser("../data/GTF/nothere.gtf", "GTF"));
}

TEST(parserGTFTest_getNextEntry, ENDL)
{
    /**< Validate opening and first value*/
    uParser ourParser(ENDL, "GTF");
    while (ourParser.eof()==false)
    {
        EXPECT_NO_THROW(ourParser.getNextEntry());
    }
}

TEST(parserGTFTest_getNextEntry, UCSCHEADER)
{
    /**< Validate opening and first value*/
    uParser ourParser(HEADER, "GTF");
    while (ourParser.eof()==false)
    {
        EXPECT_NO_THROW(ourParser.getNextEntry());
    }
}

TEST(parserGTFTest_getNextEntry, INVALIDDATA)
{
//    const int ITEMCOUNT=7;
    /**< Validate opening and first value*/
    uParser ourParser("../data/GTF/invalidsample.gtf", "GTF");
    EXPECT_ANY_THROW(ourParser.getNextEntry());
}

TEST(parserGTFTest_getNextEntry, VALIDMULTIPLE)
{
//    const int ITEMCOUNT=7;
    /**< Validate opening and first value*/
    uParser ourParser(UCSCKNOWNGENE, "GTF");
    EXPECT_NO_THROW( ourParser.getNextEntry());
    EXPECT_NO_THROW(    while(ourParser.eof()==false){
            ourParser.getNextEntry();

    });

}

TEST(parserGTFTest_getNextEntry, PARSEDVALID)
{
//    const int ITEMCOUNT=7;
    /**< Validate opening and first value*/
    uParser ourParser(UCSCKNOWNGENE, "GTF");
    EXPECT_NO_THROW( ourParser.getNextEntry());
    EXPECT_NO_THROW(    while(ourParser.eof()==false){
            ourParser.getNextEntry();

    });

}




