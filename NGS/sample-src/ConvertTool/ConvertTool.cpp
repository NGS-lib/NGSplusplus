#include <iostream>
#include <string>
#include "NGS++.h"

using namespace std;
using namespace NGS;

void usage();

int main(int argc, char* argv[])
{
    // Check arguments/
    if (argc != 5)
    {
	usage();
	return 0;
    }
    string fileName = argv[1];
    string inputType = argv[2];
    string outputName = argv[3];
    string outputType = argv[4];

    // Do the actual conversion.
    try 
    {
        uParser parser(fileName, inputType);
        uWriter writer(outputName, outputType);
        while (!parser.eof())
        {
	    writer.writeToken(parser.getNextEntry());
        }
    }
    catch(std::runtime_error& e)
    {
        cout << e.what() << endl;
    }
    catch(ugene_exception_base& e)
    {
        cout << fetchStringError(e) << endl;
    }

    return 0;
}

void usage()
{
    cout << "Usage:" << endl;
    cout << "    ConvertTool <fileName> <inputType> <outputName> <outputType>" << endl;
    cout << "    fileName: The name of the file to convert." << endl;
    cout << "    inputType: The type of the input file (i.e.: BED, SAM, WIG, etc...)." << endl;
    cout << "    outputName: The name of the output file." << endl;
    cout << "    outputType: The type of the output file." << endl;
}
