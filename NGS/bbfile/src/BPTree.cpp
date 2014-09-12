// ***************************************************************************
//   BPTree.cpp (c) 2014
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

#include "BPTree.h"
#include <stdexcept>
#include "BPTreeLeafNodeItem.h"
#include "BPTreeLeafNode.h"
#include "BPTreeChildNode.h"
#include "BPTreeChildNodeItem.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "endian_helper.h"
BPTree::BPTree()
{
    //ctor
}

BPTree::~BPTree()
{
    //dtor
}

 /*
   *    Constructor for reading in a B+ tree from a BBFile/input stream.
   *
   *    Parameters:
   *        fis - file input stream handle
   *        fileOffset - file offset to the B+ tree header
   *        isLowToHigh - indicates byte order is low to high, else is high to low
   * */
     BPTree::BPTree(std::ifstream& fis, uint64_t fileOffset) {

        // Save the seekable file handle and B+ Tree file offset
        // Note: the offset is the B+ Tree Header Table E file location
        this->fis_ = &fis;
        treeOffset_ = fileOffset;

        // read in B+ tree header - verify the B+ tree info exits
        treeHeader_ = new BPTreeHeader(*this->fis_, treeOffset_);

        // log error if header not found and throw exception
        if(!treeHeader_->isHeaderOK()){
           // int badMagic = treeHeader_->getMagic();
            throw std::runtime_error("Error reading B+ tree header: bad magic = BPTree constructor");
        }

        // assign B+ tree specifications from the header
        blockSize_ = treeHeader_->getBlockSize();
        keySize_ =  treeHeader_->getKeySize();
        valueSize_ = treeHeader_->getValSize();
        itemCount_ = treeHeader_->getItemCount();

        // populate the tree - read in the nodes
        uint64_t nodeOffset = treeOffset_ + treeHeader_->BPTREE_HEADER_SIZE;
        BPTreeNode* parentNode = NULL;  // parent node of the root is itself, or null

        // get the root node - which recursively populates the remaining nodes
        rootNode_ =  readBPTreeNode(*fis_, nodeOffset, parentNode);

    }

    /*
    *   Method returns the file input stream handle
    * */
    std::ifstream* BPTree::getFis() {
        return fis_;
    }

    /*
    *   Method returns the B+ tree file location
    * */
    uint64_t BPTree::getBPTreeOffset() {
        return treeOffset_;
    }

    /*
    *   Method returns the B+ tree header (Table E).
    * */
     BPTreeHeader* BPTree::getTreeHeader(){
        return treeHeader_;
    }

    /*
    *   Method returns the node block size (B+ order).
    * */
    uint32_t BPTree::getBlockSize() {
        return blockSize_;
    }

    /*
    *   Method returns the chromosome name key size, which is
    *   the number of valid characters for chromosome name.
    * */
     uint32_t BPTree::getKeySize() {
        return keySize_;
    }

    /*
    *   Method returns the indexing value size (currently 8).
    * */
    uint32_t BPTree::getValueSize() {
          return valueSize_;
    }

    /*
    *   Method returns the number of chromosome/contig names.
    * */
    uint64_t BPTree::getItemCount() {
        return itemCount_;
    }

    /*
    *   Method returns the number of nodes in the B+ tree.
    * */
    uint64_t BPTree::getNodeCount() {
        return nodeCount_;
    }

    /*
    *   Method returns the root node, from which all other nodes
    *   can be extracted.
    *
    *   Returns:
    *       Root node
    * */
    BPTreeNode* BPTree::getRootNode() {
        return rootNode_;
    }

    std::string BPTree::getChromosomeKey(std::string chromosome) {

        std::string key;
        if (chromosomeKeyCache_.count(chromosome)==0)
            {
            //first keySize_ chars of chrChars
            if (chromosome.size()<=keySize_)
                key =chromosome;
            else{
                key=  chromosome.substr(0,keySize_);
            }
            chromosomeKeyCache_[chromosome]=key;
        }else
        {
            key=chromosomeKeyCache_[chromosome];
        }
        return key;
    }

    /*
    *   Returns a chromosome ID  which  can be used to search for a
    *   corresponding data section in the R+ tree for data.
    *
       Parameters:
    *       chromKey - chromosome name of valid key size.
    *
    *
    *   Note: A chromosomeID of -1 means chromosome name not included in B+ tree.
    *
    * */
    uint32_t BPTree::getChromosomeID(std::string chromKey) {
         int chromosomeID;

        // Search the B+ tree to extract the Chromosome ID.
        BPTreeNode* thisNode = rootNode_;

        chromosomeID = findChromosomeID(thisNode, chromKey);

        return chromosomeID;
    }

    /*
    *   Returns a chromosome name which is the B+ key for returning the
    *   chromosome ID for lookup in the R+ tree for data.
    *
    *   Parameters:
    *       chromID - chromosome ID expected in B+ tree
    *
    *   Returns:
    *       Chromosome name key; a null string means chromosome ID not found.
    *
    * */
    std::string BPTree::getChromosomeName(uint32_t chromID) {
         std::string chromKey;

        // Search the B+ tree to extract the Chromosome ID.
        BPTreeNode* thisNode = rootNode_;

        chromKey = findChromosomeName(thisNode, chromID);

        return chromKey;
    }

    /*
    *   Method returns all chromosome key names in B+ tree.
    *
    *   Returns:
    *   Collection of all (chromosome ID, chromosome name)entries
    * */
    std::vector<std::string> BPTree::getChromosomeNames(){

        // Search the B+ tree to extract the chromosome ID.
        BPTreeNode* thisNode = rootNode_;

        std::vector<std::string> chromosomeList;

        findAllChromosomeNames(thisNode, chromosomeList);

        return chromosomeList;
    }

     /*
    *   Method returns all chromosome name, chromosome ID pairs for a given ID range.
    *
    *   Parameters:
    *       startChromID - starting ID for chromosome range expected in B+ tree
    *       endChromID - ending ID for chromosome range expected in B+ tree
    *
    *   Returns:
    *       Collection of (chromosome ID, chromosome name key) hash items;
    *       where an empty collection means ID range was not found.
    *
    * */
       std::map<uint32_t, std::string>  BPTree::getChromosomeIDMap(uint32_t startChromID, uint32_t endChromID){

        // Search the B+ tree to extract the chromosome ID.
        BPTreeNode* thisNode_ = rootNode_;

        std::map<uint32_t, std::string> chromosomeIDMap;

        findChromosomeMap(thisNode_, startChromID, endChromID, chromosomeIDMap);

        #ifdef VERBOSE
        std::cerr<< "Map count is" << chromosomeIDMap.size()<<"\n";
        #endif

        return chromosomeIDMap;
    }

    /*
    *   Method finds and returns the chromosome ID for the specified chromosome key.
    *
    *   Note: This method recursively calls itself, traversing the full B+ tree until
    *       either the chromosome name key is found and returns a valid chromosome ID,
    *       or exits with a -1 value.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       chromKey - chromosome name key of valid key size.
    *
    *   Returns:
    *       Valid chromosome ID if >= 0; else -1 for not found.
    * */
    uint32_t BPTree::findChromosomeID( BPTreeNode* thisNode, std::string chromKey){
        int chromID = -1;    // until found

        // search down the tree recursively starting with the root node
        if(thisNode->isLeaf())
        {
           int nLeaves = thisNode->getItemCount();
           for(int index = 0; index < nLeaves; ++index){
               BPTreeLeafNodeItem* leaf = (BPTreeLeafNodeItem*)thisNode->getItem(index);
               if(leaf == NULL){
                    throw std::runtime_error("Error reading B+ tree leaf nodes, corruption suspected");
               }

               // test chromosome key match
               if(leaf->chromKeysMatch(chromKey)){
                   chromID = leaf->getChromID();
                   break;
               }
               // else check next leaf
           }
        }
        else {
           // check all child nodes
           int nNodes = thisNode->getItemCount();
           for(int index = 0; index < nNodes; ++index){

               BPTreeChildNodeItem* childItem = (BPTreeChildNodeItem*)thisNode->getItem(index);
               BPTreeNode* childNode =  childItem->getChildNode();

               // check if key is in the node range
               std::string lowestKey = childNode->getLowestChromKey();
               std::string highestKey = childNode->getHighestChromKey();

               // test name key against key range
               if(chromKey!=lowestKey)  {

                    // keep going until leaf items are checked
                    chromID = findChromosomeID(childNode, chromKey);

                    // check for chromKey match
                    if(chromID >= 0)
                        break;
               }
           }
        }

        return chromID;
    }

    /*
    *   Method finds and returns the chromosome name for the specified chromosome ID.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       chromID - B+ tree chromosome ID supplied for the chromosome key
    *
    *   Returns:
    *       chromosome name if found; else a null string.
    * */
    std::string BPTree::findChromosomeName( BPTreeNode* thisNode, uint32_t chromID){

        std::string chromKey = ""; // mark unfound condition as an empty string

        // search down the tree recursively starting with the root node
        if(thisNode->isLeaf())
        {
           int nLeaves = thisNode->getItemCount();
           for(int index = 0; index < nLeaves; ++index){
               BPTreeLeafNodeItem* leaf = (BPTreeLeafNodeItem*)thisNode->getItem(index);

               if(leaf->getChromID() == chromID){ // mChromosome key match
                   chromKey = leaf->getChromKey();
                   break;
               }
               // else check next leaf
           }
        }
        else {
           // check all child nodes
           uint32_t nNodes = thisNode->getItemCount();
           for(uint32_t index = 0; index < nNodes; ++index){

               BPTreeChildNodeItem* childItem = (BPTreeChildNodeItem*)thisNode->getItem(index);
               BPTreeNode* childNode =  childItem->getChildNode();

               // check if key is in the node range
               uint32_t lowestID = childNode->getLowestChromID();
               uint32_t highestID = childNode->getHighestChromID();

               // test chromosome ID against node ID range
               if(chromID >= lowestID && chromID <= highestID) {

                    // keep going until leaf items are checked
                    chromKey = findChromosomeName(childNode, chromID);

                    // check for chromosome ID match
                    if(chromKey != "")
                        break;
               }
           }
        }

        return chromKey;
    }

    /*
    *   Method finds and returns all chromosome names in the B+ tree.
    *
    *   Note: This method calls itself recursively until the full B+ tree is traversed.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       chromosomeList - list of all chromosome names found.
    *
    *   Returns:
    *       Chromosome names found are added to the chromosome list passed in.
    * */
    void BPTree::findAllChromosomeNames( BPTreeNode* thisNode, std::vector<std::string>& chromosomeList){

        // search down the tree recursively starting with the root node
        if(thisNode->isLeaf())
        {
           // add all leaf names
           int nLeaves = thisNode->getItemCount();
           for(int index = 0; index < nLeaves; ++index){

               BPTreeLeafNodeItem* leaf = (BPTreeLeafNodeItem*)thisNode->getItem(index);
               chromosomeList.push_back(leaf->getChromKey());
           }
        }
        else {
           // get all child nodes
           int nNodes = thisNode->getItemCount();
           for(int index = 0; index < nNodes; ++index){

               BPTreeChildNodeItem* childItem = (BPTreeChildNodeItem* )thisNode->getItem(index);
               BPTreeNode* childNode = childItem->getChildNode();

               // keep going until leaf items are extracted
               findAllChromosomeNames(childNode, chromosomeList);
           }
        }
    }

    /*
    *   Method finds and returns (chromosome ID, chromosome key name) pairs for the specified ID range.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       startChromID - starting chromosome ID for the chromosome range
    *       endChromID - ending chromosome ID for the chromosome range
    *
    *   Returns:
    *       (chromosome ID, chromosome key name) items are added to the collection passed in.
    * */
    void BPTree::findChromosomeMap( BPTreeNode* thisNode, uint32_t startChromID, uint32_t endChromID,
                                        std::map<uint32_t, std::string>& chromosomeMap){
        uint32_t chromID;
        uint32_t lowestID;
        uint32_t highestID;

        // check if node is disjoint
        lowestID = thisNode->getLowestChromID();
        if(lowestID > endChromID)
            return;

        highestID = thisNode->getHighestChromID();
        if(highestID < startChromID)
            return;

        // search down the tree recursively starting with the root node
        if(thisNode->isLeaf())
        {
           int nLeaves = thisNode->getItemCount();
           for(int index = 0; index < nLeaves; ++index){

               BPTreeLeafNodeItem* leaf = (BPTreeLeafNodeItem*)thisNode->getItem(index);
               chromID = leaf->getChromID();

               // check for chromosome range match
               if( chromID >= startChromID && chromID <= endChromID ){
                    std::string key =leaf->getChromKey();
                   chromosomeMap[chromID]=key;
               }
               // leaf ID's are in ascending order; check for going out of range
              // else if(chromID > endChromID)
              //     break;
           }
        }
        else {
           // check all child nodes
           int nNodes = thisNode->getItemCount();
           for(int index = 0; index < nNodes; ++index){

               BPTreeChildNodeItem* childItem = (BPTreeChildNodeItem* )thisNode->getItem(index);
               BPTreeNode* childNode =  childItem->getChildNode();

               // check if keys are in the node range
               lowestID = childNode->getLowestChromID();
               highestID = childNode->getHighestChromID();

               // test for chromosome range intersections
               if( lowestID <= endChromID && highestID >= startChromID )
                    findChromosomeMap(childNode, startChromID, endChromID, chromosomeMap);

               // test node ID range which is always in ascending order going out of range
              // else if(lowestID > endChromID)
              //     break;   //
           }
        }
    }

    /*
    *   Method reads in the B+ tree nodes from the file, recursively.
    *
    *   Parameters:
    *       fis - file input stream handle
    *       fileOffset - file offset for B+ tree header
    *       keySize - chromosome name key size in characters
    *       parent - parent node
    *       isLowToHigh - if true, indicates byte order is low to high; else is high to low
    *
    *   Returns:
     *      Boolean which indicates if the B+ tree header was read correctly, with
    *       true for success, false for failure to find the header information.
    * */
    BPTreeNode* BPTree::readBPTreeNode(std::ifstream& fis, uint64_t fileOffset,
                                      BPTreeNode* parent){

        typedef char byte;
       // LittleEndianInputStream lbdis = null;     // low to high byte reader
       // DataInputStream bdis = null;        // high to low byte reader

        // set up for node format
     //   byte[] buffer = new byte[BPTREE_NODE_FORMAT_SIZE];
        BPTreeNode* thisNode = NULL;
        BPTreeNode* childNode = NULL;

        byte type;
        byte bval;
        uint32_t itemCount;
        uint32_t itemSize;
        bool isLeaf;

        try {

           // Read node format into a buffer
            fis.clear();
            fis.seekg(fileOffset, std::ios_base::beg);
            if (fis.eof())
                std::cerr<<"Hit end of file in seekg in readBPtreeNode\n";


            fis.read( reinterpret_cast<char*>(&type) , sizeof(byte) );
            type=endian::LittleByte(type);

           // create the B+ tree node
           if(type == 1) {
               isLeaf = true;
               thisNode = new BPTreeLeafNode(++nodeCount_);
           }
           else {
               isLeaf = false;
               thisNode = new BPTreeChildNode(++nodeCount_);
           }

            fis.read( reinterpret_cast<char*>(&bval) , sizeof(byte) );
            bval=endian::LittleByte(bval);

            fis.read( reinterpret_cast<char*>(&itemCount) , sizeof(uint16_t) );
            itemCount=endian::LittleShort(itemCount);

            // Note: B+ tree node item size is the same for leaf and child items
            itemSize =  BPTREE_NODE_ITEM_SIZE + this->keySize_;
          //  int totalSize = itemSize * itemCount;
            //byte[] itemBuffer = new byte[totalSize];

            // get the node items - leaves or child nodes
            for(uint32_t item = 0; item < itemCount; ++item) {

               // always extract the key from the node format
               std::string keychars;  // + 1 for 0 byte
               keychars.resize(keySize_);
               uint32_t index;
               for(index = 0; index < keySize_; ++index) {

                    fis.read( reinterpret_cast<char*>(&bval) , sizeof(byte) );
                    bval=endian::LittleByte(bval);

                    keychars[index] = (char)bval;
               }

                std::string key = keychars;
                std::string::iterator starti=key.begin();
                std::string::iterator endi=key.end();

                std::string::iterator end_pos = std::remove(starti,endi,' ');
                key.erase(end_pos, key.end());

               uint32_t chromID;
               uint32_t chromSize;
               uint64_t childOffset;

               if(isLeaf) {
                    fis.read( reinterpret_cast<char*>(&chromID) , sizeof(uint32_t) );
                    chromID=endian::LittleLong(chromID);

                    fis.read( reinterpret_cast<char*>(&chromSize) , sizeof(uint32_t) );
                    chromSize=endian::LittleLong(chromSize);

                    // insert leaf items
                    BPTreeLeafNodeItem* leafItem = new BPTreeLeafNodeItem(++leafCount_, key, chromID, chromSize);
                    thisNode->insertItem(leafItem);
               }
               else {
                   // get the child node pointed to in the node item
                   fis.read( reinterpret_cast<char*>(&childOffset) , sizeof(uint64_t) );
                   childOffset=endian::LittleDouble(itemCount);

                   childNode = readBPTreeNode(*fis_, childOffset, thisNode);
                   // insert child node item
                   BPTreeChildNodeItem* childItem = new BPTreeChildNodeItem(item, key, childNode);
                   thisNode->insertItem(childItem);
                }

                 fileOffset += itemSize;
           }

        }catch(...) {
           throw std::runtime_error("Error reading B+ tree node \n ");
        }

        // success: return node
        return thisNode;
   }
