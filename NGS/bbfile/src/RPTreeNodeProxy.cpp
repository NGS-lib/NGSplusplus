// ***************************************************************************
//   RPTReeNodeProxy.cpp (c) 2014
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

#include "RPTreeNodeProxy.h"

#include <stdexcept>


    RPTreeNodeProxy::~RPTreeNodeProxy()
    {
        //dtor
    }

    RPTreeNodeProxy::RPTreeNodeProxy(std::ifstream& fis, uint64_t fileOffset, int32_t chromId) {
        fis_ = &fis;
        fileOffset_ = fileOffset;
        chromId_ = chromId;
    }

    bool RPTreeNodeProxy::isLeaf() {
        throw new std::runtime_error("Not implemented -- this should never be called on a node proxy");
    }

    RPChromosomeRegion* RPTreeNodeProxy::getChromosomeBounds() {
        return NULL;
    }

    int32_t RPTreeNodeProxy::compareRegions(RPChromosomeRegion* chromosomeRegion) {
        return 0;
    }

    int32_t RPTreeNodeProxy::getItemCount() {
        return 0;
    }

    RPTreeNodeItem* RPTreeNodeProxy::getItem(int32_t index) {
        return NULL;
    }

    bool RPTreeNodeProxy::insertItem(RPTreeNodeItem* item) {
        return false;
    }

    bool RPTreeNodeProxy::deleteItem(int32_t index) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }


