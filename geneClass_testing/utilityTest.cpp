#include "utility/utility.h"
#include "gtest.h"
#include "gnuplot_i.hpp"
using namespace std;

TEST(utility, MeanTest)
{
    const vector<float> ourValues= {22, 11, 30, 0 ,50};

    EXPECT_FLOAT_EQ(utility::getMean(ourValues),22.6);

}

TEST(utility, SDTest)
{
    const vector<float> ourValues= {22, 11, 30, 0 ,50};
    EXPECT_NEAR(utility::getSd(ourValues, utility::getMean(ourValues)),19.047, 0.001);
}


TEST(utility, QuartTest)
{
    const vector<float> ourValues= {22, 11, 30, 0 ,50};


    auto quart = utility::quartilesofVector(ourValues);

    EXPECT_EQ(quart.at(0), 11);
    EXPECT_EQ(quart.at(1), 22);
    EXPECT_EQ(quart.at(2), 30);
}


TEST(utility, querySamTest)
{
    const int NONE=0;
    const int PAIRED=1;
    const int UNMAPPED=4;
    const int MINUS_STRAND=16;

   EXPECT_TRUE( utility::querySamFlag(PAIRED,SamQuery::IS_PAIRED) );
   EXPECT_FALSE( utility::querySamFlag(PAIRED,SamQuery::NEXT_UNMAPPED) );

   EXPECT_FALSE( utility::querySamFlag(MINUS_STRAND,SamQuery::FAIL_QUAL) );
   EXPECT_TRUE( utility::querySamFlag(MINUS_STRAND,SamQuery::SEQ_REV_STRAND) );

}

TEST(clustering, hausdorff)
{

    const int IDENTICAL=0;

    vector<float> vecA= {1.5,2.5,2.5,1.5};
    vector<float> vecB= {8,1,1,8};


    EXPECT_EQ(clustering::hausdorffTwoRegions(vecA,vecA), IDENTICAL);
    EXPECT_NE(clustering::hausdorffTwoRegions(vecA,vecB), IDENTICAL);


}

