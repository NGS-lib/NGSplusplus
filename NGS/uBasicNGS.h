// ***************************************************************************
// uBasicNGS.h (c) 2013
// Alexei Nordell-Markovits : Sherbrooke University
// Charles Joly Beauparlant : Laval University
//
//       This file is part of the NGS++ library.
//
//    The NGS++ library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program (lgpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// ***************************************************************************

#ifndef UBASICNGS_H_INCLUDED
#define UBASICNGS_H_INCLUDED
#include "uFormatBase.h"
#include "uFormatChrom.h"
#include "uFormatExperiment.h"
#include <memory>
namespace NGS
{

class uTags;
class uRegion;
class uGene;
class uBasicNGS: public uGenericNGS<uBasicNGS>
{
public:

    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand=StrandDir::FORWARD);
    uBasicNGS():uGenericNGS(){}
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, StrandDir pStrand, float pScore );
    uBasicNGS(std::string pChr, long long int pStart, long long int pEnd, float pScore );

    uBasicNGS(uToken pToken);
    uBasicNGS(uTags pTags);
    uBasicNGS(uRegion pRegion);
    uBasicNGS(uGene pGene);

    ~uBasicNGS(){};

    bool isEqual(const  uBasicNGS & pCompared) const;

    uBasicNGS getCopy()const;

};

class uBasicNGSChrom: public uGenericNGSChrom<uBasicNGSChrom,uBasicNGS>
{
public:
    uBasicNGSChrom():uGenericNGSChrom()
    {
    }

    uBasicNGSChrom(std::string chr):uGenericNGSChrom(chr)
    {
    }

    uBasicNGSChrom(std::string chr, long long int length):uGenericNGSChrom(chr,length)
    {
    }
  //  uBasicNGSChrom(const uGenericNGSChrom<uBasicNGSChrom,uBasicNGS> & copyCop);
    uBasicNGSChrom(const uBasicNGSChrom& initFrom);
    uBasicNGSChrom& operator=(const uBasicNGSChrom& copFrom);
    uBasicNGSChrom(const std::vector<uBasicNGS> & copyVec):uGenericNGSChrom(copyVec){};

    uBasicNGSChrom getCopy()const;
    ~uBasicNGSChrom(){};

};

class uBasicNGSExperiment : public uGenericNGSExperiment<uBasicNGSExperiment,uBasicNGSChrom,uBasicNGS>
{
    public:
    ~uBasicNGSExperiment(){};
    uBasicNGSExperiment& getCopy();
};

} // End of namespace NGS
#endif // UBASICNGS_H_INCLUDED
