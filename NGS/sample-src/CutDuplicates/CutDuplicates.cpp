#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;

int main(int argc, char **argv)
{

        int count;
        auto functOp = [&](int a)
        {

                a=a+1;
                count++;
        };


    /**< Simple parameter validation. We heavily recommand you use a library to manage your command line parameters */
    if (argc!=2)
    {
        cerr << "Please follow exactly the following signature"<<endl<<"CutDuplicates <File to eliminate duplicates>"<<endl;
        return 0;
    }
    ifstream firstStream;
    {
        string firstPath=argv[1];
        try
        {
            /**< Load file in a corresponding stream. Will throw if path is invalid */
            utility::loadStream(firstPath,firstStream);
        }
        catch(...)
        {
            /**< If invalid filePath, error message and abort. At this point, we have not validated the file format. */
            cerr << "Error loading file. Check your path"<<endl;
            return 0;
        }
    }
    /**< Catch block, for any number of possible errors. This will also catch Parser errors */
  try {
     uBasicNGSExperiment loadedFile;
     /**< Load every item from our first file.
     This will throw if the data is poorly formated
     */
     vector<string> customColumns={"CHR","START_POS","END_POS"};
    uParser customParser(&firstStream,customColumns);

    loadedFile.loadWithParser(customParser);

    /**< Sort every element by region */
    loadedFile.sortSites();
    /**< Create our writer. Note that this means we are not keeping the original lines.
    We pass cout so that it will write to standard output and can be pipped
    */
    uWriter bedWriter(&cout,"BED6");
    /**< We declare a simple lambda function we will be using to check for overlaps
    Note that [&] in the signature means will be be capturing by reference whatever outside
    object we care to use.
    */
    set<pair<long long, long long>> countSet;
    auto functElimDuplicate=[&](uBasicNGSChrom & chrom)
    {
    countSet.clear();
     /**< Already sorted */
    for (auto it = chrom.begin(); it!=chrom.end(); it++)
        {
            if (!countSet.count(pair<long long,long long>(it->getStart(),it->getEnd())))
            {
                countSet.insert(pair<long long,long long>(it->getStart(),it->getEnd()));
                it->writeToOutput(bedWriter);
            }
        }
    };

    /**< We passe our previously defined lambda to this function
        It will loaded the data from our second BED file and item by item
        will run the passed lambda function. Note this means our second BED file is never
        loaded into memory.
     */
        loadedFile.applyOnAllChroms(functElimDuplicate);


    }
    /**< Global exception handling. Practicalloy, this should probably be split so as to
         manage parser errors seperately.
     */
    catch(ugene_exception_base & e)
    {
        cout << fetchStringError(e)<<endl;
    }

}
