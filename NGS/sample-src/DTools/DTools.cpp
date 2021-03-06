#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include "NGS++.h"

using namespace std;
using namespace NGS;


/** \brief This program load a sam file then outputs density information, either in wig format or in "BIN" format.
           It uses NGS++ mostly as a datastructure to query the data with examples of function binding, using the Parser
           and the getSubsetCount function.

           Writing could have been done with the Writer function, but is kept manual for now.
*/


struct wigData
{
    std::string chr;
    long long int position;
    long long int span;
    float value;
};

vector<wigData> vectorToWig(vector<long int> densityVector, string chrom);
void writeWig(const vector<wigData> & ourData, std::ostream& out);
void densityFromFile(int argc, char* argv[]);
void writeBinDensity( uTagsChrom& tagChrom, std::ostream& out, int binSize);
void writeDensityFromFile(const uTagsChrom& tagChrom, std::ostream& out);

int main(int argc, char* argv[])
{
   string firstArg="";
    densityFromFile(argc, argv);
    return 0;
}

void densityFromFile(int argc, char* argv[])
{
try{
    string samPath="";
    string outputPath="";
    int binSize =0;

    /**< Parsing input parameters, we highly recommend you use TCLAP instead ( see the extra section of the NGS++ library) */
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i],"-s")==0)
        {
            samPath = argv[i + 1];
        }
        else if (strcmp(argv[i],"-o")==0)
        {
            outputPath= argv[i + 1];
        }
        else if (strcmp(argv[i],"-b")==0)
        {
            binSize= std::atoi(argv[i + 1]);
        }
    }

    /**< If a mandatory parameter is missing, fail */
    if ((samPath.size()==0))
    {
        cerr<<"Program signature is  -s <SamPath>  -o [OutputPath] -b [binSize = Optional]\n";
        return;
    }

    /**< If path to file is invalid, fail */
    ifstream inputStream(samPath);
    if( !( inputStream))
    {
        cerr <<"Invalid file path for Sam File, failling"<<endl;
        return;
    }
    /**< We declare the standard unit for SAM type files */
    uTagsExperiment tagExp;
    /**< Load our data from Stream */
    uParser samParser(samPath,"SAM");
    tagExp.loadWithParser(samParser);
    /**< Query header info */
    try {
    auto chrVector=samParser.getHeaderParamVector(header_param::CHR);
    auto chrSizeVector=samParser.getHeaderParamVector(header_param::CHR_SIZE);

    if ((int)chrVector.size()==0){
        cerr <<"Missing scaffold header information, please validate your SAM header is complete ( or present )"<<std::endl;
        return;
    }

    for(size_t i=0; i<chrVector.size();i++)
    {
        if (tagExp.isChrom(chrVector.at(i)))
            tagExp.setChrSize(chrVector.at(i),std::stoi(chrSizeVector.at(i)));
    }
    }
    catch(...){
        cerr <<"Failed while loading header information,please validate your SAM header";
        return;
    }

    /**< Sort our data per start bp position */
    tagExp.sortSites();

    std::ofstream outputOS;
    if (outputPath.size()!=0)
    {
        outputOS.open(outputPath);
    }
    /**< If we did not have an output path, output to standard output */
    std::ostream & outFile = ((outputPath.size()!=0) ? outputOS : std::cout);
    if (binSize)
    {
    /**<  We create a function object that takes a uTagsChrom as only parameter*/
    /**< Note that we must wrap our stream in a reference wrapper or it will modify a copy ( and as copying a stream is disabled by default, will not compile */
       std::function<void(uTagsChrom& tagChrom)> writeBind =std::bind(&writeBinDensity,std::placeholders::_1,ref(outFile),binSize);
       /**< Apply this function on all chromosomes, this will write the density of every bin. */
        tagExp.applyOnAllChroms(writeBind);
    }
    else
    {
        for (auto chromIT= tagExp.begin(); chromIT!=tagExp.end(); chromIT++)
        {
           /**< Idem as above, but for wig rather then binned output */
         std::function<void(const uTagsChrom&)> writeDensityWig=std::bind(&writeDensityFromFile,std::placeholders::_1,ref(outFile));
         tagExp.applyOnAllChroms(writeDensityWig);
        }
        }
    }
    /**< Various exception management */
    catch (elem_throw & e)
    {
        const uTags* errorTagPoint;
        if( ( errorTagPoint=boost::get_error_info<tag_error>(e) ) )
        {
            cerr <<" We crashed working on this Sam Tag" <<endl;
        }
        if (std::string const * ste =boost::get_error_info<string_error>(e) )
            {
            cerr << "Trace of crash" <<endl;
            cerr << *ste;
        }
    }
    catch(ugene_exception_base &e)
    {
        cerr << "Caught exception from Parser object, failling" <<endl;
        cerr <<fetchStringError(e)<<endl;
    }
    catch(std::exception &e)
    {
        cerr << "Caught standard exception, failling" <<endl;
        cerr <<e.what()<<endl;
    }
}

/**< For every bin, get the count of element who's START overlaps this bin and write it out */
void writeBinDensity( uTagsChrom& tagChrom, std::ostream& out, int binSize)
{
    vector<long int> densityValues;
    densityValues.resize(tagChrom.getChromSize());
    string chrName= tagChrom.getChr();

    for (int j=0; j <((int)densityValues.size()/binSize); j++ )
    {
        int start=j*binSize;
        int end=j*binSize+binSize;
        int tagcount=0;
        /**< This returns the count of elements who's start is between start and end bp */
        /**< See manual for how this could have been used for different type of subsets */
        tagcount= tagChrom.getSubsetCount(start,end);
        /**< Write the bind */
        if (tagcount)
            out<<chrName << "\t" <<start << "\t" <<end <<"\t" << "." << "\t" <<tagcount <<endl;
    }
}

void writeDensityFromFile(const uTagsChrom& tagChrom, std::ostream& out)
{
    vector<long int> densityValues;
    /**< Create a density map for the entire chromosome. Using a sparse Vector would be more efficient */
    densityValues.resize(tagChrom.getChromSize());
    try {
        /**< Query every site to get data */
        const_cast<uTagsChrom*> (&tagChrom)->applyOnAllSites([&] (uTags Elem)
        {
            for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
            {
                if (i<(int)densityValues.size())
                    densityValues.at(i)++;
                else{
                    cerr << "Skipping the following tags as over chromSize of " <<densityValues.size() <<endl;
                    Elem.print(cout);
                }
            }
        }
                                                        );
        auto ourWig=vectorToWig(densityValues,tagChrom.getChr() );
        writeWig(ourWig, out);
    }
    catch(std::exception &e)
    {
        cerr << "Catching and re-throwing from writeDensityFromFile" <<endl;
        cerr << "Failed on chr "<<tagChrom.getChr() <<endl;
        cerr << e.what()<<endl;
        throw e;
    }

}
/**< Using our Parser class would be even simple, kept as demonstration purpose */
void writeWig(const vector<wigData> & ourData, std::ostream& out)
{
try {
    const string step="variableStep chrom=";
    const string span=" span=";
    auto curSpan= 0;
    for(auto& wigValue:ourData )
    {
        if (wigValue.span!=curSpan)
        {
            out <<step<<wigValue.chr<<span<<wigValue.span<<endl;
            curSpan=wigValue.span;
        }
        out << wigValue.position << "\t" << wigValue.value << endl;
    }
}
catch(std::exception &e)
{
    cerr << "Catching and re-throwing from writeWig" <<endl;
}
}

vector<wigData> vectorToWig(vector<long int> densityVector, string chrom)
{
    try {
    wigData tempData;
    tempData.position=0;
    tempData.chr=chrom;
    tempData.value=densityVector.at(0);
    vector<wigData> returnVec;
    long long int span=0;
    long long int pos=0;;
    for(auto& densityValue:densityVector )
    {
        if (densityValue!=tempData.value)
        {
            tempData.span=span;
            returnVec.push_back(tempData);
            /**< New data */

            tempData.value = densityValue;
            tempData.position=pos;
            span=1;
        }
        else
            span++;

        /**< Next position */
        pos++;
    }
    if (tempData.span>1)
        returnVec.push_back(tempData);

    return returnVec;
    }
    catch(std::exception &e)
    {
        cerr <<"Throwing from vectorToWig"<<endl;
        throw e;

    }
}
