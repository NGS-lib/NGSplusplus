#ifndef UGENE_H_INCLUDED
#define UGENE_H_INCLUDED

#include "uFormats.h"
#include "chromConstants.h"
namespace NGS {

//This is a strongly typed class that represents a gene entry.
//Our problem.. is that a gene is many, many things.
//Let us hard code our mandatory entriesà that we take from the GenePred format.
//But we want this to be extensible, for future analyse.
//http://www.sequenceontology.org/gff3.shtml
class uGeneFeature::uGenericNGS {

  //Easier to keep
  struct exon{
      int exonStart;
      int exonEnd;
  };
  private:

  //
  std::string ID;
  //Strand is one of these questions we should ask...should this be lower down?
  char strand;

//Coding start, coding stop.
  int cdsStart;
  int cdsEnd;
/*
  long int txStart;
  long int txEnd;


  long int cdsStart;
  long int cdsEnd;
  int exonCount;
  int FiveUtr;
  int ThreeUtr;
  std::vector<exon> exonList;
  std::string proteinID;
  std::string alignID;
 */
//Fonctions
 // std::string getChromfromint(int chromNum);
 // int getChromfromString(std::string);

  public:

  uGene();
  ~uGene();



//  friend std::ostream &operator<<(std::ostream &out, const uGene& ourGene); // output

  void writeToOuput(std::ostream &out);

  long int distanceFromTss( int pos);
  bool isInGene(int pos);
  bool overlapsGene(int start, int end);
  bool isInExon(int pos);
  bool overlapsTss(int start, int end);
  bool isTss(int pos, int tssSize);
  bool isProximalProm(int pos, int proximalSize);


  //Accesors
  std::string getChrom(){return chrom;};
  int getTssStart(){return txStart;};
  //int getChromNum(){return chromNum;};

  void loadfromUCSCString(std::string geneInfo);

};


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
} // End of namespace NGS

