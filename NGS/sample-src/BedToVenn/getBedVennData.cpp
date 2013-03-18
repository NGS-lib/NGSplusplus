#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;

/**< Input is sorted */

uBasicNGSExperiment returnOverlapParts(const uBasicNGSExperiment& set1, const uBasicNGSExperiment& set2)
{
    /**< For every chrom */
    uBasicNGSExperiment returnExp;
    for(auto chromItr1= set1.begin();chromItr1!=set1.end(); chromItr1++ )
    {
        /**< If other set has chrom */
        if (set2.isChrom(chromItr1->first))
        {
            /**< Compare data in chrom */
            auto pChrComp = set2.getpChrom(chromItr1->first);
            const uBasicNGSChrom * pChrCur= &(chromItr1->second);
            /**< Compare every item. Brute force */
            {
                /**< First chrom loop */
                for (auto curChromItr= pChrCur->begin();curChromItr!=pChrCur->end();curChromItr++ )
                {
                    /**< Compare to every other element */
                    for (auto compChromItr= pChrComp->begin();compChromItr!=pChrComp->end();compChromItr++ ){
                        if (compChromItr->getStart()>curChromItr->getEnd())
                            break;

                    if(utility::isOverlap(curChromItr->getStart(),curChromItr->getEnd(), compChromItr->getStart(),compChromItr->getEnd()))
                        returnExp.addData(curChromItr->returnOverlapping(*compChromItr));

                    }
                }
            }
        }
    }
    return returnExp;
}

void generateAndOutput(const pair<std::string,uBasicNGSExperiment> & createdSet,  vector<pair<std::string,uBasicNGSExperiment>> ::iterator curItr,vector<pair<std::string,uBasicNGSExperiment>> ::iterator endItr  ,uWriter& pWriter)
{
    uBasicNGSExperiment compareResult =returnOverlapParts(createdSet.second, curItr->second);
/**< If no elements, stop comparing, pointless */
    if(compareResult.count())
    {
          std::string newID= createdSet.first +"_"+curItr->first;
          cout << ("Compare "+newID)<<endl;
          compareResult.writeWithWriter(pWriter);
       /**< Call all possible combinations */
       curItr++;
       while(curItr!=endItr){
            generateAndOutput(make_pair(newID,compareResult),curItr,endItr,pWriter);
            curItr++;
       }
    }
}


int main(int argc, char **argv)
{
    /**< Simple parameter validation. We heavily recommand you use a library to manage your command line parameters */
    if (argc<=2)
    {
        cerr << "Please follow exactly the following signature"<<endl<<"<File1 to compare> <File2 to Compare>"<<endl;
        return 0;
    }
    vector<pair<std::string,uBasicNGSExperiment> > vecBedExp;
    vecBedExp.resize(argc-1);
    try
    {
        for (int i=1; i<argc; i++ )
        {
            {
             string path=argv[i];
             ifstream curStream;
             utility::loadStream(path,curStream);
             uParser curParser(&curStream,"BED");
             vecBedExp.at(i-1).first=path;
             vecBedExp.at(i-1).second.loadWithParser(curParser);
             vecBedExp.at(i-1).second.sortSites();
            }
        }
    }
    catch(...)
    {
        /**< If invalid filePath, error message and abort. At this point, we have not validated the file format. */
        cerr << "Error loading files. Check your path"<<endl;
        return 0;
    }


    /**< Catch block, for any number of possible errors. This will also catch Parser errors */
  try {

    /**< Compare every example with each other. */

        uWriter coutWriter(&cout, "BED4");

        auto secundItr=(vecBedExp.begin()+1);
        auto endItr= (vecBedExp.end());

         while(secundItr!=endItr){
            generateAndOutput(vecBedExp.at(0),secundItr,endItr,coutWriter);
            secundItr++;
       }


//        for (size_t i=0; i<vecBedExp.size(); i++)
//        {
//            for(size_t y=i+1; y<vecBedExp.size(); y++)
//            {
//                auto compareResult=returnOverlapParts(vecBedExp.at(i).second,vecBedExp.at(y).second);
//                if (compareResult.count()){
//                        std::string newID= vecBedExp.at(i).first +"_"+vecBedExp.at(y).first;
//                    cout << ("Compare "+newID)<<endl;
//                    compareResult.writeWithWriter(coutWriter);
//                    vecBedExp.push_back()
//                }
//            }
//        }

     }
    catch(ugene_exception_base & e)
    {
        cout << fetchStringError(e)<<endl;
    }

}
