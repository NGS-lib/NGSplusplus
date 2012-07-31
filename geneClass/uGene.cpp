#include "uGene.h"
#include "limits.h"

uGeneRef::uGeneRef() {  // default constructor


}

uGeneRef::~uGeneRef() {  // default destructor

}

//Read from a standard UCSC Known Gene DB
void uGeneRef::loadfromUCSCDB(ifstream& geneStream){

string geneInfo;
uGene tempGene;


//Check if header
    getline(geneStream,geneInfo);

    if (geneInfo[0]!='#'){
          geneStream.seekg(0);
    }


 while(!geneStream.eof()){


     getline(geneStream, geneInfo);
    //geneStream.getline(geneInfo,500);

    tempGene.loadfromUCSCString(geneInfo);

    addgene(tempGene);

    if (geneStream.eof())
        break;
 }

}

//Add a Gene to our data
  void uGeneRef::addgene(uGene ourgene){


    UgeneVec* pVecGene;

    pVecGene=&(mapOfGenes[ourgene.getChrom()]);


   pVecGene->push_back(ourgene);

  }


//Output in UCSC tab delimited format
  void uGeneRef::outputChromtoFile(string chrom, ofstream& out ){

    UgeneVec* pVecGene;


     pVecGene=&(mapOfGenes[chrom]);

    for (unsigned int i=0 ;i< pVecGene->size(); i++){
        pVecGene->at(i).writeToOuput(out);
    }
  }


 void uGeneRef::outputAlltoFile(ofstream& out){

    //Typedef used
    GeneVectorMap::iterator iterMap;

     for (iterMap = mapOfGenes.begin(); iterMap != mapOfGenes.end(); ++iterMap) {
        outputChromtoFile(iterMap->first, out);
    }
 }

//Our data has to be sorted.
uGene uGeneRef::findNearestGene(string chrom, int pos)
{
    //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return pGene;

    //GotoReplace by limits.h
    long int closest=2147483647;
    long int current=0;


    for (unsigned int i=0; i< pVecGene->size(); i++ ){

        current=pos-pVecGene->at(i).getTssStart();
            if (abs(current) < abs(closest)){
                closest = current;
                pGene= pVecGene->at(i);
            }
    }


return pGene;
}

//Our data has to be sorted.
uGene uGeneRef::findNearestPriorGene(string chrom, int pos){

     //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return pGene;

    long int closest=2147483647;
    long int current=0;


    for (unsigned int i=0; i< pVecGene->size(); i++ ){
        current=pos-pVecGene->at(i).getTssStart();
            if  ((abs(current) < abs(closest))&&(current>=0)){
                closest = current;
                pGene= pVecGene->at(i);
            }
    }

return pGene;
}

//Our data has to be sorted.
uGene uGeneRef::findNearestAfterGene(string chrom, int pos){

   //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return pGene;

    long int closest=2147483647;
    long int current=0;

    for (unsigned int i=0; i< pVecGene->size(); i++ ){
        current=pos-pVecGene->at(i).getTssStart();
            if ( (abs(current) < abs(closest))&&(current<=0)  ){
                closest = current;
                pGene= pVecGene->at(i);
            }
    }

return pGene;

}

vector<uGene> uGeneRef::findAllRegion(string chrom,int pos, int size){

    vector<uGene> GeneVector;
    int current;

  //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return GeneVector;

    for (unsigned int i=0; i< pVecGene->size(); i++ ){
        current=pos-pVecGene->at(i).getTssStart();
            if ( (abs(current) < size)&&(current<=0)  ){
                pGene= pVecGene->at(i);
               GeneVector.push_back(pGene);

            }
    }


    return GeneVector;
}

 bool uGeneRef::isTss(std::string chrom,int pos, int size){
    //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return 0;



    for (unsigned int i=0; i< pVecGene->size(); i++ ){
       if (pVecGene->at(i).isTss(pos, size)){
            return true;
       }
    }

 return false;

 }


 bool uGeneRef::isProximalProm(std::string chrom,int pos, int size){

   //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return 0;



    for (unsigned int i=0; i< pVecGene->size(); i++ ){
       if (pVecGene->at(i).isProximalProm(pos, size)){
            return true;

       }
    }

 return false;
 }

   bool uGeneRef::isGene(std::string chrom,int pos){

   //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return 0;

    for (unsigned int i=0; i< pVecGene->size(); i++ ){
       if (pVecGene->at(i).isInGene(pos)){
            return true;
       }
    }
 return false;
}

 bool uGeneRef::overlapsTss(std::string chrom,int start, int end){

  //Get the Vectors of chrom
  UgeneVec* pVecGene;

  uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return 0;


    for (unsigned int i=0; i< pVecGene->size(); i++ ){
       if (pVecGene->at(i).overlapsTss(start, end)){
            return true;
       }
    }

 return false;

 }

//Do we overlap a gene
 bool uGeneRef::overlapsGene(std::string chrom, int start, int end){

 //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return 0;

    for (unsigned int i=0; i< pVecGene->size(); i++ ){
       if (pVecGene->at(i).overlapsGene(start,end)){
            return true;
       }
    }

return false;
 }

//Is pos in an Exon
bool  uGeneRef::isExon(std::string chrom,int pos){

 //Get the Vectors of chrom
      UgeneVec* pVecGene;

      uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return 0;

    for (unsigned int i=0; i< pVecGene->size(); i++ ){
       if (pVecGene->at(i).isInExon(pos)){
            return true;
       }
    }
 return false;



}

//Distance of our position to nearest TSS
long int uGeneRef::distanceFromClosestTss(std::string chrom,int pos){

 UgeneVec* pVecGene;
 int currentvalue=LONG_MAX;
  uGene pGene;

    pVecGene=&(mapOfGenes[chrom]);
    if (pVecGene->size()==0)
        return 0;


    for (unsigned int i=0; i< pVecGene->size(); i++ ){
       if (abs((pVecGene->at(i).distanceFromTss(pos)))<currentvalue){
            currentvalue=(pVecGene->at(i).distanceFromTss(pos));
       }
    }

 return currentvalue;

}


using namespace std;

uGene::uGene() {  // default constructor


}

uGene::~uGene() {  // default destructor

}

//Load info from a standard UCSC string
void uGene::loadfromUCSCString(string geneInfo){

stringstream geneInfostream;
char tempString[255];

    geneInfostream.str(geneInfo);


    geneInfostream>>chrom;
    geneInfostream>>GeneName;


   // chromNum=getChromfromString(chrom);
    geneInfostream>>strand;
    geneInfostream>>txStart;
    geneInfostream>>txEnd;
    geneInfostream>>cdsStart;
    geneInfostream>>cdsEnd;
    geneInfostream>>exonCount;
    exonList.resize(exonCount);
    for (int i=0 ; i <exonCount; i++){
        geneInfostream.getline(tempString,255,',');
        exonList.at(i).exonStart= atoi(tempString); }
    for (int i=0 ;i <exonCount; i++){
        geneInfostream.getline(tempString,255,',');
        exonList.at(i).exonEnd= atoi(tempString); }

    //Discard next tab
    getline(geneInfostream, proteinID,'\t');
    //Get until next
    getline(geneInfostream, proteinID,'\t');
    //Idem
    getline(geneInfostream, alignID,'\n');
 //   geneInfostream>>proteinID;
  //  geneInfostream>>alignID;

}





 void uGene::writeToOuput(std::ostream &out)
{

    out << chrom << "\t" << GeneName << "\t" << strand << "\t" << txStart << "\t" <<
    txEnd << "\t" <<  cdsStart << "\t" <<  cdsEnd << "\t" <<  exonCount<<"\t";

    for (int i=0 ; i <exonCount; i++){
        out  << exonList.at(i).exonStart<<",";
    }
    out <<"\t";
    for (int i=0 ;i <exonCount; i++){
         out << exonList.at(i).exonEnd<< "," ;
    }
    out<< "\t" << proteinID <<"\t"<< alignID <<"\n";
    return;
}

//This assumes the correct chromosome has been checked.
  long int uGene::distanceFromTss(int pos){

    return (pos-txStart);

  };


  bool uGene::isInGene(int pos){

    if (utility::isInInterval(pos, txStart, txEnd))
        return true;

    return false;
  }



  bool uGene::isInExon(int pos){

   for (int i=0 ; i <exonCount; i++){

       if (utility::isInInterval(pos,exonList.at(i).exonStart,exonList.at(i).exonEnd)){
            return true;
       }

    }

return false;
}

//Is the point in a pre-determined size around our TSS
bool uGene::isTss(int pos, int tssSize){

    int txStartSpace, txEndSpace;

    txStartSpace= txStart- tssSize;
        if (txStartSpace<0)
            txStartSpace=0;
    txEndSpace=( txStart+ tssSize);


     if (utility::isInInterval(pos,txStartSpace,txEndSpace)){
            return true;
       }

   return false;
  }


//WE consider the promoter to end at the transcription start site.
bool uGene::isProximalProm(int pos, int proximalSize){

  int PromStart,PromEnd;
    PromStart = txStart;
    PromEnd =txStart;

    PromStart= txStart- proximalSize;
        if (PromStart<0)
            PromStart=0;

    if (utility::isInInterval(pos, PromStart,PromEnd)){
            return true;
       }

   return false;

}

//Does our site overlap our Gene
 bool uGene::overlapsGene(int start, int end){


    //We start in Gene
    if ((start> txStart)&&(start < txEnd))
        return true;

    //We end in Gene
    if ((end> txStart)&&(end < txEnd))
        return true;

    //We envelop the gene
    if ((start< txStart)&&(end > txEnd))
        return true;

return false;
 }

bool uGene::overlapsTss(int start, int end){


    if ((start <= txStart) && ( end >= txStart ))
        return true;

return false;
}

