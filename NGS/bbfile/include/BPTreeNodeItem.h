#ifndef BPTREENODEITEM_H_INCLUDED
#define BPTREENODEITEM_H_INCLUDED

// ***************************************************************************
//   BPTreeNodeItem.h (c) 2014
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
#include <string>
/*
*   BPTreeNodeItem interface for storage of B+ tree node item information.
*
*   Note: The alpha-numeric key string is used for positional insertion of
*    node items and searching of the B+ tree.
* */

class BPTreeNodeItem{


    public:
    virtual ~BPTreeNodeItem() {}

    // Returns the child node item or leaf item index in the B+ tree.
    virtual uint64_t getItemIndex()=0;

    // Identifies the item as a leaf item or a child node item.
    virtual bool isLeafItem()=0;

    // Returns key used to position the item in parent node item list.
    virtual std::string getChromKey()=0;

    // Returns true if keys match, returns false if keys do not match.
    virtual bool chromKeysMatch(std::string chromKey)=0;

};

#endif // BPTREENODEITEM_H_INCLUDED
