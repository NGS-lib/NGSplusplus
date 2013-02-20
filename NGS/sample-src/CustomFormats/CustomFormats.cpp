#include <iostream>
#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;
int main(int argc, char **argv)
{
     /**< Simple parameter validation. We heavily recommand you use a library to manage your command line parameters */

    ifstream firstStream;
    ifstream secondStream;
    {
        string firstPath=argv[1];
        /**< Load each file in a corresponding stream. Will throw if path is invalid */
        utility::loadStream(firstPath,firstStream);
        utility::loadStream(firstPath,secondStream);

    }
    /**< This determines the order that we will store the read data */
    std::vector<std::string> columnFormat={"START_POS","END_POS","SCORE","CHR"};
    /**< As we are passing a vector rather then a type, the parser automatically sets itself to "CUSTOM" */
    /**< As an optional third parameter, we pass ',', this is our custome delimiter. By default, \t is used */
    uParser customParser(&firstStream, columnFormat,',');

    uWriter bedWriter(&cout,"BED6");
    /**< Write every read line as Bed 6 format */
    while (!customParser.eof()){
        bedWriter.writeToken(customParser.getNextEntry());
    }

    cout <<"Same data in custom Format" <<endl;
    /**< Load a second parser.  */
    uParser secondParser(&secondStream, columnFormat,',');

    std::vector<std::string> outputFormat={"SCORE","START_POS","END_POS","CHR"};

    uWriter cusTomWriter(&cout,outputFormat,'-');
    while (!(secondParser.eof())){
        cusTomWriter.writeToken(secondParser.getNextEntry());
    }

}

