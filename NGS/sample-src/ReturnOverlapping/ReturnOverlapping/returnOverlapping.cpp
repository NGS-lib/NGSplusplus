#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;

int main(int argc, char **argv)
{
    /**< Simple parameter validation. We heavily recommand you use a library to manage your command line parameters */
    if (argc!=3)
    {
        cerr << "Please follow exactly the following signature"<<endl<<"returnOverlapping <first File> <second File>"<<endl;
        return 0;
    }


    ifstream firstStream;
    ifstream secondStream;
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
     uBasicNGSExperiment loadedFile;
     /**< Load every item from our first file.
     This will throw if the data is poorly formated
     */
    loadedFile.loadWithParser(firstStream,"BED");
    /**< Create our writer. Note that this means we are not keeping the original lines.
    We pass cout so that it will write to standard output and can be pipped
    */
    uWriter bedWriter(&cout,"BED4");

    /**< We declare a simple lambda function we will be using to check for overlaps
    Note that [&] in the signature means will be be capturing by reference whatever outside
    object we care to use.
    */
    auto functOverlap=[&](uBasicNGS & item)
    {
        /**< Validate the item we are comparing is from a scaffold present in our previously loaded EXP */
        if (loadedFile.isChrom(item.getChr()))
        {
            /**< Created a dummy Chrom structure containing only our item */
            /**< There are more efficient way's, but this is simple */
            uBasicNGSChrom compareChrom(item.getChr());
            compareChrom.addData(item);

            /**< Get the pointer to the Chrom scaffold from our loaded EXP */
            auto pChrom= loadedFile.getpChrom(item.getChr());
            /**< If our item overlaps at least one item from our EXP, we write it to output. */
            if (pChrom->getOverlappingCount(compareChrom)>0)
                item.writeToOutput(bedWriter);
        }
    };
    /**< We passe our previously define lambda to this function
        It will loaded the data from our second BED file and item by item
        will run the passed lambda function. Note this means our second BED file is never
        loaded into memory.
     */
    loadedFile.loadWithParserAndRun(secondStream,"BED",functOverlap);

    }
    /**< Global exception handling. Practicalloy, this should probably be split so as to
         manage parser errors seperately.
     */
    catch(ugene_exception_base & e)
    {
        cout << fetchStringError(e)<<endl;
    }

}
