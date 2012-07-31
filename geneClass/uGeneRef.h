#ifndef UGENEREF_H_INCLUDED
#define UGENEREF_H_INCLUDED

#include "uGene.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>


#include <map>

 using namespace std;
class uGeneRef{


  private:


typedef std::vector<uGene> UgeneVec;
typedef map<string, vector<uGene> >  GeneVectorMap;

 std::string assemblyname;
 std::vector<uGene> geneList;

//Map on the chromosome string
 std::map<string, vector<uGene> > mapOfGenes;
 //Utility file

  public:

  uGeneRef();
  ~uGeneRef();

  void addgene(uGene ourgene);
  int count(){return geneList.size();};

  uGene findNearestGene(string chrom, int pos);
  uGene findNearestPriorGene(string chrom, int pos);
  uGene findNearestAfterGene(string chrom, int pos);
  bool  isTss(std::string chrom,int pos, int size);
  bool  isProximalProm(std::string chrom,int pos, int size);
  bool  isGene(std::string chrom,int pos);
  bool  overlapsTss(std::string chrom,int start, int end);
  bool  overlapsGene(std::string chrom, int start, int end);
  bool  isExon(std::string chrom,int pos);
  long int distanceFromClosestTss(std::string chrom,int pos);

  std::vector<uGene> findAllRegion(string chrom,int pos, int size);

 void outputChromtoFile(string chrom, ofstream& );
 void outputAlltoFile(ofstream&);


  //std::vector<uGene> findRegion(string chrom,int pos, int size);
  void loadfromUCSCDB(ifstream& geneStream);


};



#endif // UGENEREF_H_INCLUDED
