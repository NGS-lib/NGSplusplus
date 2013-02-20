#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;


void uBasicToCout(const uBasicNGS & item){
    std::cout <<"Peak Called:"<<"\t"<< item.getChr() <<"\t"<< item.getStart() <<"\t"<< item.getEnd() <<std::endl;
}

void uGeneToCout(const uGene & item){
      std::cout << item.getChr() <<"\t"<< item.getStart() <<"\t"<< item.getEnd() <<"\t"<<item.getID() <<"\t"<<item.getTranscript()<<std::endl;
}

int main(int argc, char **argv)
{
    /**< Simple parameter validation. We heavily recommand you use a library to manage your command line parameters */
    if (argc!=4)
    {
        cerr << "Please follow exactly the following signature"<<endl<<"geneOverlap <MACS bed file> <UCSC genePred reference> <File to Write Overlaps>"<<endl;
        return 0;
    }
    ifstream firstStream;
    ifstream secondStream;
    string thirdPath=argv[3];
    {
        string firstPath=argv[1];
        string secondPath = argv[2];

        try
        {
            /**< Load each file in a corresponding stream. Will throw if path is invalid */
            utility::loadStream(firstPath,firstStream);
            utility::loadStream(secondPath,secondStream);
        }
        catch(...)
        {
            /**< If invalid filePath, error message and abort. At this point, we have not validated the file format. */
            cerr << "Error loading files. Check your path"<<endl;
            return 0;
        }
    }
    /**< Catch block, for any number of possible errors. This will also catch Parser errors */
  try {
     uGeneExperiment UCSCRefData;
     /**< Load every item from our UCSC ref. This will automatically parse in the EXONS and coding region
     This will throw if the data is poorly formated
     */
    UCSCRefData.loadWithParser(secondStream,"GENEPRED");
    UCSCRefData.sortSites();
    string thirdPath=argv[3];
    uWriter bedWriter(thirdPath, "BED6");
    /**< We declare a lambda function we will be using to check for overlaps and write the info
    */
    auto functOverlap=[&](uBasicNGS item)
    {
        /**< Validate the item we are comparing is from a scaffold present in our previously loaded EXP */
        if (UCSCRefData.isChrom(item.getChr()))
        {

            /**< Get the pointer to the Chrom scaffold from our loaded EXP */
            auto pChrom= UCSCRefData.getpChrom(item.getChr());
            /**< If our item overlaps at least one item from our EXP, we write it to output. */
            if (pChrom->getSubsetCount(item.getStart(),item.getEnd())>0)
            {
                /**< Write the item, since it does overlap */
                bedWriter.writeToken(item.createToken());
                uBasicToCout(item);
                auto overlappingItems= pChrom->getSubset(item.getStart(),item.getEnd());
                cout << "The above called region overlaps the following "+to_string(overlappingItems.count())+" genes and the following features of each gene" <<endl;
                for (auto itr=overlappingItems.begin();itr!=overlappingItems.end();itr++)
                {
                    uGeneToCout(*itr);
                   // std::cout <<"Feature count is "<<itr->featureCount()<<std::endl;
                    for (auto featureItr= itr->featureBegin(); featureItr!=itr->featureEnd();featureItr++)
                    {
                        if (utility::isOverlap(featureItr->getStart(),featureItr->getEnd(),item.getStart(),item.getEnd()))
                            std::cout << item.getChr() <<"\t"<< featureItr->getStart() <<"\t"<< featureItr->getEnd()<<"\t"<<featureStr(featureItr->getType())<<"\t"<<std::endl;
                       // else
                       //     std::cout<< "Not overlapping"<<item.getChr() <<"\t"<< featureItr->getStart() <<"\t"<< featureItr->getEnd()<<"\t"<<featureStr(featureItr->getType())<<"\t"<<std::endl;
                    }
                }
                cout <<std::endl;
                /**< For each gene body overlapped, write */


                 item.writeToOutput(bedWriter);
            }

        }
    };
    /**< We passe our previously define lambda to this function
        It will loaded the data from our second BED file and item by item
        will run the passed lambda function. Note this means our second BED file is never
        loaded into memory.
     */
        uBasicNGSExperiment interExp;
        vector<std::string> format{"CHR","START_POS","END_POS","JUNK","SCORE"};
        uParser customeParse(&firstStream,format);
        interExp.loadWithParser(customeParse);
        interExp.applyOnSites(functOverlap);
    }
    /**< Global exception handling. Practicalloy, this should probably be split so as to
         manage parser errors seperately.
     */
    catch(ugene_exception_base & e)
    {
        cout << fetchStringError(e)<<endl;
    }

}
