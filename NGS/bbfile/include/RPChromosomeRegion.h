#ifndef RPCHROMOSOMEREGION_H
#define RPCHROMOSOMEREGION_H



// ***************************************************************************
//   RPChromosomeRegion.h (c) 2014
//   Copyright @ Alexei Nordell-Markovits : Sherbrooke University
//
//    This file is part of the BWReader library.
//
//    The BWReader library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU  General Public License
//    along with this program (gpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// ***************************************************************************

// The code structure and logic is based on the original IGV BBFileReader.
// The original code can be found here:
//https://github.com/broadinstitute/IGV
// The code was originally released under the LGPL 2.1
//http://www.opensource.org/licenses/lgpl-2.1.php).
// Our thanks to the IGV team for making the code available.

#include "stdint.h"

class RPChromosomeRegion
{
    public:
        RPChromosomeRegion();
        virtual ~RPChromosomeRegion();
        RPChromosomeRegion(uint32_t startChromID, uint32_t startBase, uint32_t endChromID, uint32_t endBase);
        RPChromosomeRegion(const RPChromosomeRegion* region);

        uint32_t getStartChromID() { return startChromID_; }
        void setStartChromID(uint32_t val) { startChromID_ = val; }
        uint32_t getStartBase() { return startBase_; }
        void getStartBase_(uint32_t val) { startBase_ = val; }
        uint32_t getEndChromID() { return endChromID_; }
        void setEndChromID(uint32_t val) { endChromID_ = val; }
        uint32_t getEndBase() { return endBase_; }
        void setEndBase(uint32_t val) { endBase_ = val; }

        uint32_t compareRegions(RPChromosomeRegion* testRegion);
        bool equals(RPChromosomeRegion* testRegion);


        bool containedIn(RPChromosomeRegion* testRegion);
        bool intersectsBelow(RPChromosomeRegion* testRegion);
        bool intersectsAbove(RPChromosomeRegion* testRegion);
        bool disjointBelow(RPChromosomeRegion* testRegion);
        bool disjointAbove(RPChromosomeRegion* testRegion);
        RPChromosomeRegion* getExtremes(RPChromosomeRegion* testRegion);
        void expand(RPChromosomeRegion* testRegion);

    protected:

        uint32_t startChromID_;
        uint32_t startBase_;
        uint32_t endChromID_;
        uint32_t endBase_;
};

#endif // RPCHROMOSOMEREGION_H
