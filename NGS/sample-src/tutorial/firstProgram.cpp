#include "NGS++.h"
using namespace std;
using namespace NGS;

int main(int argc, char* argv[]) {
	string pathToSam("someValidPath"); // Obviously this path is not valid, for educationnal purpose only
	uTagsExperiment ourExperiment, newExp;
	uParser samParser(pathToSam, "SAM");
	ourExperiment.loadWithParser(samParser);
	auto chromSiteOver200 = [](const uTagsChrom chr) {return ( chr.count()>200) ;};
	map<string,uTagsChrom> chrMap = ourExperiment.getSpecificChroms(chromSiteOver200);
	uTagsExperiment newExperiment;
	for (auto x: chrMap) {
		newExperiment.addData(x.second);
	}

	uWriter bedWriter(&cout,"BED4");
	ourExperiment.writeWithWriter(bedWriter);
}
