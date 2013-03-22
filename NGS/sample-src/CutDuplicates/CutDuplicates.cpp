#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;

int main(int argc, char **argv)
{
    /**< Simple parameter validation. We heavily recommand you use a library to manage your command line parameters */
    if (argc!=3)
    {
        cerr << "Please follow exactly the following signature"<<endl<<"<File to eliminate duplicates> [BED/GFF/SAM]"<<endl;
        return 0;
    }
    string fileType=argv[2];
    if (fileType=="GFF")
        fileType="UCSCGFF";

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
     uTagsExperiment loadedFile;
     /**< Load every item from our first file.
     This will throw if the data is poorly formated
     */
    uParser fileParser(&firstStream,fileType);

    /**< Create our writer. Note that this means we are not keeping the original lines.
    We pass cout so that it will write to standard output and can be pipped
    */
    if (fileType=="BED")
       uWriter bedWriter(&cout,"BED6");
    else
         uWriter bedWriter(&cout,fileType);





    map<std::string,map<long int, set<long int> >> duplicateMap;

    long int removedCount=0;
    while(fileParser.eof()==false){
        uBasicNGS myTags(fileParser.getNextEntry());
        if (duplicateMap[myTags.getChr()][myTags.getStart()].count(myTags.getEnd())==0)
        {
            duplicateMap[myTags.getChr()][myTags.getStart()].insert(myTags.getEnd());
            cout <<fileParser.getPreviousRaw()<<"\n";
        }else
        {
            removedCount++;
        }
    }

    std::cerr <<"Removed "<<removedCount<<std::endl;
    }
    /**< Global exception handling. Practicalloy, this should probably be split so as to
         manage parser errors seperately.
     */
    catch(ugene_exception_base & e)
    {
        cerr << fetchStringError(e)<<endl;
    }

}

/**<THIS IS ANOTHER WAY OF DOING THE TASk. KEPT FOR PURPOSE OF EXAMPLE  */


    /**< We declare a simple lambda function we will be using to check for overlaps
    Note that [&] in the signature means will be be capturing by reference whatever outside
    object we care to use.
    */

//    loadedFile.loadWithParser(fileParser);
//
//    /**< Sort every element by region */
//    loadedFile.sortSites();

//    long long removedCount=0;
//    map<std::string,map<pair<long,long>,std::string>> complexMap;
//    set<pair<long long, long long>> countSet;
//    auto functElimDuplicate=[&](uTagsChrom & chrom)
//    {
//    countSet.clear();
//     /**< Already sorted */
//    for (auto it = chrom.begin(); it!=chrom.end(); it++)
//        {
//            if (!countSet.count(pair<long long,long long>(it->getStart(),it->getEnd())))
//            {
//                countSet.insert(pair<long long,long long>(it->getStart(),it->getEnd()));
//                it->writeToOutput(bedWriter);
//            }
//            else{
//                removedCount++;
//            }
//        }
//    };
//
//    /**< We passe our previously defined lambda to this function
//        It will loaded the data from our second BED file and item by item
//        will run the passed lambda function. Note this means our second BED file is never
//        loaded into memory.
//     */
//        loadedFile.applyOnAllChroms(functElimDuplicate);
