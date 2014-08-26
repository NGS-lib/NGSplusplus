#ifndef BPTREE_H
#define BPTREE_H


// ***************************************************************************
//   BPTree.h (c) 2014
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


#include "BPTreeNode.h"
#include "BPTreeHeader.h"
#include "stdint.h"
#include <string>
#include <map>
#include <vector>
class BPTree
{
    public:
        BPTree();
        virtual ~BPTree();

        static const  uint32_t BPTREE_NODE_FORMAT_SIZE = 4;   // node format size
        static const  uint32_t BPTREE_NODE_ITEM_SIZE = 8;     // Plus keySize to be added

        protected:
        private:
        // B+ tree access variables   - for reading in B+ tree nodes from a file
        std::ifstream* fis_;      // file handle - BBFile input stream
        uint64_t treeOffset_;         // mChromosome B+ tree file offset
        BPTreeHeader* treeHeader_; // B+ tree header (Table E for BBFile)

        // B+ tree organizational variables  - derived from Table E
        uint32_t blockSize_;     // number of children per block
        uint32_t keySize_;       // character size of primary key
        uint32_t valueSize_;     // number of bytes in value being indexed
        uint64_t itemCount_;    //  number of contig/mChromosome items in tree

        // B+ tree nodal variables
        BPTreeNode* rootNode_;   // B+ tree root node
        uint64_t nodeCount_;        // number of nodes defined in the B+ tree
        uint64_t leafCount_;        // number of leaves in the B+ tree


        public:

         BPTree(std::ifstream& fis, uint64_t fileOffset);
         std::ifstream* getFis() ;
         uint64_t getBPTreeOffset();

        std::map<std::string, std::string> chromosomeKeyCache_;
        BPTreeHeader* getTreeHeader();
        uint32_t getBlockSize() ;
        uint32_t getKeySize();
        uint32_t getValueSize();
        uint64_t getItemCount();
        uint64_t getNodeCount();
        uint32_t getChromosomeID(std::string chromKey);
        std::string getChromosomeName(uint32_t chromID);
        std::vector<std::string> getChromosomeNames();
        std::map<uint32_t, std::string>  getChromosomeIDMap(uint32_t startChromID, uint32_t endChromID);
        uint32_t findChromosomeID( BPTreeNode* thisNode, std::string chromKey);
        std::string findChromosomeName( BPTreeNode* thisNode, uint32_t chromID);
        void findAllChromosomeNames( BPTreeNode* thisNode, std::vector<std::string>& chromosomeList);
        void findChromosomeMap( BPTreeNode* thisNode, uint32_t startChromID, uint32_t endChromID,
                                        std::map<uint32_t, std::string>& chromosomeMap);

        BPTreeNode* readBPTreeNode(std::ifstream& fis, uint64_t fileOffset, BPTreeNode* parent);
        BPTreeNode* getRootNode();
        std::string getChromosomeKey(std::string chromosome);

};

#endif // BPTREE_H
