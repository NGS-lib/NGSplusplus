#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;

int main(int argc, char **argv)
{
    /**< Simple parameter validation. We heavily recommand you use a library to manage your command line parameters */
    if (argc!=2)
    {
        cerr << "Please follow exactly the following signature"<<endl<<"<Input bed>"<<endl;
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
    /**< Declare our Bedgraph writer and Bed reader */
    uParser bedParser(&firstStream,"BED");

    uBasicNGSExperiment toSortData;
    toSortData.loadWithParser(bedParser);
    toSortData.sortSites();
    uWriter bedGraphWriter(&cout,"BEDGRAPH");
    bedGraphWriter.writeHeader();
    toSortData.writeWithWriter(bedGraphWriter);
  //
  //  while (bedParser.eof()==false){

  //      bedGraphWriter.writeToken(bedParser.getNextEntry());

  //  }
    /**< Write our UCSC header */
//    UCSCHeader::getUCSCBrowserLine

    /**< Write our bedgrpah */

    }
    /**< Global exception handling. Practicalloy, this should probably be split so as to
         manage parser errors seperately.
     */
    catch(ugene_exception_base & e)
    {
        cout << fetchStringError(e)<<endl;
    }

}
