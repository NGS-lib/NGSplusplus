#include "NGS++.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include <time.h>
#include "gtest.h"
#include <ctime>

time_t tstart, tend;
using namespace NGS;
using namespace std;

TEST(uTagsTest, TimeTestingRawParser){


	time_t tstart, tend;
	tstart = time(0);
	try {
	uParser samParser("../data/SAM/H2AZFinal_downsample.sam","SAM");

	while(samParser.eof()==false)
		samParser.getNextEntry();
	}
	catch(ugene_exception_base &e){
	cout << fetchStringError(e)<<endl;
	}

tend = time(0);
cout << "It took with Parser" << difftime(tend, tstart) << " second(s)." << endl;

}

TEST(uTagsTest, TimeTestingInnerParser){

	time_t tstart, tend;
	tstart = time(0);
	try {
    ifstream ourStream;
    utility::loadStream("../data/SAM/H2AZFinal_downsample.sam",ourStream);
    uTagsExperiment ourSam;
    ourSam.loadFromSamWithParser("../data/SAM/H2AZFinal_downsample.sam");
	}
	catch(ugene_exception_base &e){
	cout << fetchStringError(e)<<endl;
	}

tend = time(0);
cout << "It took with inner Regex" << difftime(tend, tstart) << " second(s)." << endl;

}

TEST(uTagsTest, TimeTestingInner){


	time_t tstart, tend;
	tstart = time(0);
	try {
    ifstream ourStream;
    utility::loadStream("../data/SAM/H2AZFinal_downsample.sam",ourStream);
    uTagsExperiment ourSam;
    ourSam.loadFromSam(ourStream);

	}
	catch(ugene_exception_base &e){
	cout << fetchStringError(e)<<endl;
	}

tend = time(0);
cout << "It took with Inner" << difftime(tend, tstart) << " second(s)." << endl;
}
