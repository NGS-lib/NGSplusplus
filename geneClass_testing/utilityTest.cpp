#include "utility.h"
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


TEST(clustering, hausdorff)
{

    const int IDENTICAL=0;

    vector<float> vecA= {1.5,2.5,2.5,1.5};
    vector<float> vecB= {8,1,1,8};


    EXPECT_EQ(clustering::hausdorffTwoRegions(vecA,vecA), IDENTICAL);
    EXPECT_NE(clustering::hausdorffTwoRegions(vecA,vecB), IDENTICAL);



    //  EXPECT_EQ(clustering::hausdorffTwoRegions(vecA,vecB), IDENTICAL);

    //  cerr << "Clustering is" << clustering::hausdorffTwoRegions(vecA,vecB) <<endl;
}



TEST(clustering, Alignement)
{

    const int IDENTICAL=0;

  //  vector<float> vecA= {0,0,0,1};
  //  vector<float> vecB= {8,8,8,8};

   // vector<float> vecTryA {0,33,33,33,26,26,26,26,26,27,24,22,25,20,21,24,24,24,24,25,23,23,23,23,23,22,20,20,19,19,15,16,14,18,18,18,18,18,18,18,17,17,17,17,17,20,20,20,19,19,19,15,15,14,11,11,11,11,9,9,9,10,10,11,11,13,13,13,13,13,12,12,8,8,8,8,8,8,8,8,8,8,9,8,5,5,6,7,7,7,8,8,8,8,8,8,9,9,10,10,11,11,10,11,10,11,11,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,15,14,13,14,17,16,16,16,16,16,16,15,16,15,15,13,13,13,13,12,11,11,9,9,9,10,12,14,14,16,17,17,17,17,17,17,16,16,17,16,16,16,15,12,16,16,16,17,17,18,18,17,17,17,17,20,20,19,19,20,20,20,20,20,19,16,15,15,13,12,12,12,12,12,12,13,13,12,12,13,13,13,13,9,9,9,8,9,8,9,9,9,9,9,6,8,8,8,7,7,7,9,15,16,16,15,15,15,16,16,16,16,16,16,18,18,19,19,18,18,18,18,18,21,21,21,20,20,19,19,19,19,19,19,17,17,19,19,19,19,17,11,10,10,10,10,10,9,9,10,11,11,11,8,8,7,7,7,7,7,7,7,4,4,4,4,4,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,9,9,9,12,12,12,11,10,12,12,12,12,12,12,12,12,13,13,15,15,15,15,15,16,13,13,13,13,13,13,13,13,11,11,11,10,11,11,11,10,10,10,7,7,7,8,8,7,7,7,7,7,7,8,8,7,7,5,5,5,5,7,6,6,6,6,6,6,6,6,6,6,6,6,7,6,6,6,6,6,8,9,9,9,8,8,7,7,7,7,7,7,6,6,7,8,8,8,8,8,6,7,7,9,9,9,11,13,13,13,13,13,13,12,12,12,12,13,13,11,11,11,11,11,11,11,11,11,11,11,11,11,11,10,9,9,9,10,10,10,9,10,8,8,9,7,6,6,6,6,6,6,6,6,6,6,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,
  //                        };
  //  vector<float> vecTryB {0,0,0,0,0,0,19,20,20,20,18,18,18,17,17,19,19,19,19,19,18,18,18,18,18,18,16,15,15,15,10,10,7,7,7,7,7,7,7,7,6,6,6,6,5,7,6,6,5,5,5,4,4,4,2,2,2,3,3,3,4,4,4,4,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,4,5,5,5,5,5,7,7,7,8,9,9,8,8,8,8,8,8,9,9,12,12,12,12,12,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,19,17,18,18,17,16,16,16,16,16,15,15,15,15,14,10,10,10,10,10,11,12,22,24,24,26,26,26,25,25,25,25,25,25,26,25,25,25,26,22,22,21,21,21,21,22,24,24,25,25,26,26,25,25,25,25,25,26,27,25,24,17,15,15,13,13,13,13,13,13,13,14,14,13,13,14,14,13,14,14,16,16,16,16,15,13,13,12,12,11,12,13,13,14,14,14,13,12,13,13,10,12,12,13,16,16,16,16,16,16,20,20,20,20,19,19,19,18,18,27,28,28,28,28,29,29,30,30,30,29,28,28,29,29,29,29,29,28,29,30,28,28,27,24,24,24,27,27,27,23,23,23,23,23,23,23,22,23,12,11,11,11,11,11,11,10,10,11,11,11,11,10,11,11,12,12,12,11,10,10,10,11,11,11,11,8,8,8,7,7,7,7,8,8,9,9,9,9,9,9,9,9,8,8,9,9,8,8,8,9,8,8,8,7,7,7,8,8,8,8,7,7,7,8,8,10,10,10,11,11,11,10,10,9,10,11,11,11,11,13,13,14,14,13,13,13,13,13,12,13,12,12,12,12,12,11,11,11,11,11,11,11,10,10,9,9,9,8,9,10,10,10,11,10,8,8,8,8,6,7,6,6,6,6,6,7,8,8,7,7,7,7,7,7,7,7,7,7,9,9,9,9,9,8,8,8,8,7,6,6,6,5,5,5,5,5,5,5,4,5,6,6,7,7,6,5,5,5,5,5,5,5,5,6,6,6,7,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,
  //                        };

    vector<float> vecTryA{ 2,5,1,5,1,1,1,1};
    vector<float> vecTryB{ 2,1,1,1,5,1,1,5};

  //  double maxA = *std::max_element(vecA.begin(), vecA.end());
  //  double maxB = *std::max_element(vecB.begin(), vecB.end());
  //  double maxAbs = std::max(maxA,maxB);

    //  EXPECT_EQ(clustering::align_distance(vecB,vecB,maxA,20), IDENTICAL);
    //  EXPECT_NE(clustering::align_distance(vecA,vecB,maxAbs,20), IDENTICAL);

   cerr << "Distance Big  is" << clustering::align_distance(vecTryA,vecTryB,0,1) <<endl;



  //  vector<float> distanceVec;
  //  for (int i=0; i<(int)vecTryA.size(); i++){
   //     distanceVec.push_back(clustering::align_distance(vecTryA,vecTryB,0,i) );
  //      std::cerr<< distanceVec.back()<<std::endl;
   //     }

  /*  Gnuplot g1("lines");
    g1.set_style("lines");
    g1.plot_x(distanceVec,"at i distance< no penalty");
     distanceVec.clear();
    for (float k=0.001; k<0.5; k+=0.01)
    {
        for (int i=2; i<(int)vecTryA.size(); i++)
            distanceVec.push_back(clustering::align_distance(vecTryA,vecTryB,0,i,k) );
       // g1.reset_all();
        g1.plot_x(distanceVec,"Penalty ="+utility::numberToString(k) );
        distanceVec.clear();
    }



    g1.showonscreen();

    utility::pause_input();
    cerr << "Distance smaller is" << clustering::align_distance(vecA,vecB,0,1) <<endl; */
}




TEST(clustering, Size)
{

    vector<float> size_test;
    size_test.resize(30000);
    for (int i=0; i< 30000; i++)
        size_test.at(i)=i;
    //  std::cerr << "Big test" <<std::endl;
    clustering::align_distance(size_test,size_test,30001,100);


}
