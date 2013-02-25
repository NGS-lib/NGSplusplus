#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;

int main(int argc, char **argv)
{
    cout << UCSCHeader::getUCSCPositionLine("chr1",100,2000);
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseTypeMult::HIDE_ALL,"I","WISH","TO","FIND");
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseTypeMult::PACK_ALL,"I","WISH","TO","PACK");
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseTypeMult::DENSE_ALL,"I","WISH","TO","DENSE");
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseTypeMult::FULL_ALL,"I","WISH","TO","FULL","MANY","FULL");
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseTypeMult::SQUISH_ALL,"I","WISH","TO","SQUISH","MANY");


    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseType::HIDE);
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseType::PACK);
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseType::DENSE);
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseType::FULL);
    cout << UCSCHeader::getUCSCLine(UCSCHeader::UCSCBrowseType::SQUISH);


    /**< Simple parameter validation. We heavily recommand you use a library to manage your command line parameters */
    if (argc!=4)
    {
        cerr << "Please follow exactly the following signature"<<endl<<"fileToSplit"<<endl;
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

    std::string type = argv[2];
    std::string strID = argv[3];
    if ((type!="BED")&&(type!="GFF")&&(type!="GTF")){
         /**< If invalid filePath, error message and abort. At this point, we have not validated the file format. */
            cerr << "Please elect type from BED/GFF/GTF"<<endl;
            return 0;
    }


    /**< Catch block, for any number of possible errors. This will also catch Parser errors */
  try {
         /**< Declare our Parser */
        uParser typeParser(&firstStream,type);
        /**< Splitting a file by line is a very simple operation. Simpler to access the raw file without passing by an EXP structure */
        while(!typeParser.eof())
        {
            uToken token =typeParser.getNextEntry();
            std::string curChr= token.getParam(token_param::CHR);
            ofstream curStream(curChr+"_"+strID,ios::out | ios::app );
            curStream <<typeParser.getPreviousRaw()<<endl;
        }
    }
    /**< Global exception handling. Practicalloy, this should probably be split so as to
         manage parser errors seperately.
     */
    catch(ugene_exception_base & e)
    {
        cout << fetchStringError(e)<<endl;
    }

}
