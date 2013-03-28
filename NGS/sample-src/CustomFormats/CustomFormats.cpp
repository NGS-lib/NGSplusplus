#include <iostream>
#include <iostream>
#include "NGS++.h"
using namespace std;
using namespace NGS;

vector<string> validChoices={"START_POS","END_POS","CHR","SCORE","STRAND","SEQ_NAME","JUNK"};
bool isStringinVector(const string & value, vector<string> vectorCheck)
{
     if ( find (vectorCheck.begin(), vectorCheck.end(),value)==vectorCheck.end())
        return false;
    else
        return true;
}

bool isValid(const string &value)
{
    return isStringinVector(value, validChoices);

}
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
    vector<string> loadedOrder;
    char strSep = (argv[2][0]);

    if (strSep=='t')
        strSep='\t';
    if (strSep=='s')
        strSep=' ';

    for (int i=3; i<argc;i++)
    {
        std::string curValue(argv[i]);

        if (isValid(curValue) && isStringinVector(curValue,loadedOrder)==false )
            loadedOrder.push_back(curValue);
    }

    /**< As we are passing a vector rather then a type, the parser automatically sets itself to "CUSTOM" */
    uParser customParser(&firstStream, loadedOrder,strSep);

    uWriter bedWriter(&cout,"BED6");
    /**< Write every read line as Bed 6 format */
    while (!customParser.eof()){
        bedWriter.writeToken(customParser.getNextEntry());
    }
}

