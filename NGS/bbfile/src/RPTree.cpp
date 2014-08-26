// ***************************************************************************
//   RpTree.cpp (c) 2014
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

#include "RPTree.h"
#include <vector>
#include <stdexcept>
#include <string>
#include "endian_helper.h"
#include <stdlib.h>
#include <sstream>


#include "RPTreeChildNode.h"
#include "RPTreeLeafNode.h"
#include "RPTreeLeafNodeItem.h"
#include "RPTreeChildNodeItem.h"
#include "RPTreeNodeProxy.h"



RPTree::RPTree()
{
    //ctor
}

RPTree::~RPTree()
{
    //dtor
}



    /*
   * Constructor for reading in a B+ tree from a BBFile/input stream.
   * */
    /*
    *   Constructor for R+ chromosome data locator tree
    *
    *   Parameters:
    *       fis - file input stream handle
    *       fileOffset - location for R+ tree header
    *       isLowToHigh - binary values are low to high if true; else high to low
    *       uncompressBuffSize - buffer size for decompression; else 0 for uncompressed data
    * */
    RPTree::RPTree(std::ifstream& fis, uint64_t fileOffset, uint32_t uncompressBuffSize, bool forceDescend) {

        // save the seekable file handle  and B+ Tree file offset
        // Note: the offset is the file position just after the B+ Tree Header
        // mBBFis = fis;
        rpTreeOffset_ = fileOffset;
        uncompressBuffSize_ = uncompressBuffSize;

        // read in R+ tree header - verify the R+ tree info exits
        rpTreeHeader_ = new RPTreeHeader(fis, rpTreeOffset_);

        // log error if header not found and throw exception
        if (!rpTreeHeader_->isHeaderOK()) {
            throw std::runtime_error("Error reading R+ tree header: bad magic = ");
        }

        // assigns R+ tree organization from the header
        order_ = rpTreeHeader_->getBlockSize();
        chromosomeBounds_ = new RPChromosomeRegion(rpTreeHeader_->getStartChromID(), rpTreeHeader_->getStartBase(),
                rpTreeHeader_->getEndChromID(), rpTreeHeader_->getEndBase());

        // populate the tree - read in the nodes
        uint64_t nodeOffset = rpTreeOffset_ + rpTreeHeader_->getHeaderSize();
      //  RPTreeNode* parentNode = NULL;      // parent node of the root is itself, or null

        // start constructing the R+ tree - get the root node
        rootNode_ = readRPTreeNode(&fis, nodeOffset, forceDescend);
    }

    /*
     *  Constructs an R+ Tree which conforms to the supplied information
     *      order -  the items per node factor, sometimes called the m factor,
     *      where any node must have at least m/2 items and no more than m items.
     *      keySize - the number of significant bytes in a item key.
     * */

    RPTree::RPTree(uint32_t order) {

        // R+ tree node specification
        this->order_ = order;
        // Note: acknowledge no bounds specified as a null  object
        chromosomeBounds_ = NULL;

    }


    /*
    *   Method finds the bounding chromosome region in R+ tree for a chromosome ID range.
    *
    *   Parameters:
    *       startChromID - start chromosome for the region
    *       endChromID - end chromosome for the region
    *
    *   Returns:
    *       Region which bounds the extremes of chromosome ID range
    * */

    RPChromosomeRegion* RPTree::getChromosomeRegion(uint32_t startChromID, uint32_t endChromID) {

        RPChromosomeRegion* region;

        // Search the R+ tree to extract the chromosome region.
        RPTreeNode* thisNode = rootNode_;
        RPChromosomeRegion* seedRegion = NULL;  // null until a chromosome match

        region = findChromosomeRegion(thisNode, startChromID, endChromID, seedRegion);

        return region;
    }

    /*
    *   Method returns list of all chromosome regions found for the chromosome ID range.
    *
    *   Returns:
    *       List of all chromosome regions in the chromosome ID range.
    * */

    std::vector<RPChromosomeRegion*> RPTree::getAllChromosomeRegions() {

        // Search the R+ tree to extract the chromosome regions
        RPTreeNode* thisNode = rootNode_;

        std::vector<RPChromosomeRegion*> regionList;

        findAllChromosomeRegions(thisNode, regionList);

        return regionList;
    }

    /*
    *   Method extracts a hit list of chromosome data file locations for a specified chromosome region.
    *
    *   Parameters:
    *       chromosomeRegion - chromosome region for feature extraction consists of:
    *           startChromID - start chromosome ID for region
    *           mStartBase - starting base for data extraction
    *           endChromID - end chromosome ID for region
    *           mEndBase - ending base for data extraction
    *       contained - if true indicates all returned data must be
    *           completely contained within the extraction region;
    *           else if false, returns all intersecting region features
    *
    *   Note: The selection region will be limited to accommodate  maxLeafHits; which terminates
    *       selection at the leaf node at which maxLeafHits is reached. Total number of selected
    *       items may exceed maxLeafHits, but only by the number of leaves in the cutoff leaf node.
    *
    *   Returns:
    *       List of chromosome leaf items which identify file locations for bed data
    *       of a chromosome region, or a sub-region subject to maxLeafHits.
    *
    *       Check returned leaf item bounds for cutoff limits on selection region due to maxLeafHits.
    * */

   std::vector<RPTreeLeafNodeItem*> RPTree::getChromosomeDataHits(RPChromosomeRegion* selectionRegion, bool contained) {

        std::vector<RPTreeLeafNodeItem*> leafHitItems;

        // check for valid selection region - return empty collection if null
        if (selectionRegion == NULL)
            return leafHitItems;

        // limit the hit list size
        /*
        if(maxLeafHits > 0)
            mMaxLeafHits = maxLeafHits;
        else
            mMaxLeafHits = mRPTreeHeader_->getBlockSize();
        */

        findChromosomeRegionItems(rootNode_, selectionRegion, leafHitItems);

        return leafHitItems;
    }


    /*
    *   Method finds and returns the bounding chromosome region for the specified
    *   chromosome ID range.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       startChromID - start chromosome ID for region
    *       endChromID - end chromosome ID for region
    *       region  - leaf region contains extremes for given chromosome ID range
    *
    *   Note: region grows recursively to match the extremes found for the
    *   specified chromosome ID range.  Starting base comes from the startChromID
    *   match and ending base comes form the endChromID match.
    *
    *   Returns:
    *       Chromosome region if found in the R+ tree node passed in;
    *       else null region
    * */

    RPChromosomeRegion* RPTree::findChromosomeRegion(RPTreeNode* thisNode, uint32_t startChromID, uint32_t endChromID, RPChromosomeRegion* region) {

//        int hitValue;
        RPChromosomeRegion* bounds;

        // search down the tree recursively starting with the root node
        if (thisNode->isLeaf()) {
            int nLeaves = thisNode->getItemCount();
            for (int index = 0; index < nLeaves; ++index) {

                RPTreeLeafNodeItem* leaf = (RPTreeLeafNodeItem*) thisNode->getItem(index);

                // get leaf region bounds
                bounds = leaf->getChromosomeBounds();

                // test this leaf's chromosome ID's for chromosome hit, then include its base bounds
                if ( ( startChromID >= bounds->getStartChromID() && startChromID <= bounds->getEndChromID() ) ||
                        ( endChromID >= bounds->getStartChromID() && endChromID <= bounds->getEndChromID()) ){

                    // Note: need a start region before comparing other regions for extremes
                    if (region == NULL)
                        region = new RPChromosomeRegion(bounds); // seed extreme region
                    else
                        region = region->getExtremes(bounds); // update seed extreme region
                }
            }
        } else {
            // check all child nodes
            uint32_t nNodes = thisNode->getItemCount();
            for (uint32_t index = 0; index < nNodes; ++index) {

                RPTreeChildNodeItem* childItem = (RPTreeChildNodeItem*) thisNode->getItem(index);

                // get bounding region and compare chromosome ID's
                bounds = childItem->getChromosomeBounds();

                // test node chromosome ID range for any leaf hits for either startChromID or endChromID
                if (  ( startChromID >= bounds->getStartChromID() && startChromID <= bounds->getEndChromID() ) ||
                       ( endChromID >= bounds->getStartChromID() && endChromID <= bounds->getEndChromID())  ){

                    RPTreeNode* childNode = childItem->getChildNode();
                    region = findChromosomeRegion(childNode, startChromID, endChromID, region);
                }

                // check next node
            }
        }

        return region;
    }

    /*
    *   Method finds and returns all chromosome regions in the R+ chromosome data tree.
    *
    *   parameters:
    *       thisNode - tree node to start search
    *       chromosomeList - list of all chromosome names found.
    *
    *   Returns:
    *       Adds chromosome regions if found in the chromosome region list passed in.
    * */

    void RPTree::findAllChromosomeRegions(RPTreeNode* thisNode,
                                          std::vector<RPChromosomeRegion* > regionList) {

        // search down the tree recursively starting with the root node
        if (thisNode->isLeaf()) {
            int nLeaves = thisNode->getItemCount();
            for (int index = 0; index < nLeaves; ++index) {

                RPTreeLeafNodeItem* leaf = (RPTreeLeafNodeItem*) thisNode->getItem(index);

                // add all leaf regions
                RPChromosomeRegion* region = leaf->getChromosomeBounds();
                regionList.push_back(region);
            }
        } else {
            // get all child nodes
            int nNodes = thisNode->getItemCount();
            for (int index = 0; index < nNodes; ++index) {

                RPTreeChildNodeItem* childItem = (RPTreeChildNodeItem*) thisNode->getItem(index);
                RPTreeNode* childNode = childItem->getChildNode();

                findAllChromosomeRegions(childNode, regionList);
            }
        }

    }

    /*
    *   Method returns an array of chromosome leaf items for the chromosome test region.
    *
    *   Note: At the leaf item level, any hit is valid. Use of contained is exercised
    *       when the leaf item data region is examined.
    *   Parameters:
    *
    *       thisNode - tree node to start search
    *       testRegion - bounding region for feature extraction consists of:
    *           start chromosome ID for region
    *           starting base for data extraction
    *           end chromosome ID for region
    *           ending base for data extraction
    *       leafHitItems - array containing previous leaf hit items
    *
    *   Note: leaf hit items will be limited to leaves in the current leaf node,
    *       once the maximum number of leaf hits mMaxLeafHits is reached.
    *
    *   Returns:
    *       ArrayList of leaf hit items containing updated hit regions and data offsets;
    *       else an empty array if hit regions not found.
    * */

    void RPTree::findChromosomeRegionItems(RPTreeNode* thisNode, RPChromosomeRegion* selectionRegion,
                                           std::vector<RPTreeLeafNodeItem*> &leafHitItems) {

        int hitValue;

        // check for valid selection region - ignore request if null
        if (selectionRegion == NULL)
            return;

        // check if node is disjoint
        hitValue = thisNode->compareRegions(selectionRegion);
        if (std::abs(hitValue) >= 2)
            return;

        // search down the tree recursively starting with the root node
        if (thisNode->isLeaf()) {
            uint32_t nLeaves = thisNode->getItemCount();
            for (uint32_t index = 0; index < nLeaves; ++index) {
                RPTreeLeafNodeItem* leafItem = (RPTreeLeafNodeItem*) thisNode->getItem(index);

                // compute the region hit value
                hitValue = leafItem->compareRegions(selectionRegion);

                // select contained or intersected leaf regions - item selection is by iterator
                 if (std::abs(hitValue) < 2){
                    leafHitItems.push_back(leafItem);
                }

                // ascending regions will continue to be disjoint so terminate nodal search
                else if (hitValue > 1)
                    break;

                // check next leaf
            }
        } else {
            // check all child nodes
            uint32_t nNodes = thisNode->getItemCount();
            for (uint32_t index = 0; index < nNodes; ++index) {
                RPTreeChildNodeItem* childItem = (RPTreeChildNodeItem* ) thisNode->getItem(index);

                // check for region intersection at the node level
                hitValue = childItem->compareRegions(selectionRegion);

                // test this node and get any leaf hits; intersections and containing
                 if (std::abs(hitValue) < 2) {
                    RPTreeNode* childNode = childItem->getChildNode();
                    findChromosomeRegionItems(childNode, selectionRegion, leafHitItems);
                }

                // ascending regions will continue to be disjoint so terminate nodal search
                else if (hitValue > 1)
                    break;
            }
        }
    }

    /*
    *   Method reads in the R+ tree nodes recursively.
    *
    *   Note: If node is a child node, the node is examined recursively,
    *       until the leaves are found.
    *
    *   Parameters:
    *       fis - file input stream handle
    *       fileOffset - file location for node specification (Table L)
    *       parent - parent node of this node
    *       isLowToHigh - indicates formatted data is low to high byte order if true;
    *           else is high to low byte order
    *
    *   Returns:
    *       A tree node, for success, or null for failure to find the node information.

    * */

    RPTreeNode* RPTree::readRPTreeNode(std::ifstream* fis, uint64_t fileOffset, bool forceDescend) {

        RPTreeNode* thisNode = NULL;
        typedef int8_t byte;
        try {

            // Read node format into a buffer
            fis->clear();
            fis->seekg(fileOffset);

            // find node type
            byte type;

            fis->read( reinterpret_cast<char*>(&type) , sizeof(byte) );
            type=endian::LittleByte(type);

            bool isLeaf;
            uint32_t itemSize;


            //TODO Check the creator of RPTreeNodes, make sure parent value is called and whatbnot
            if (type == 1) {
                isLeaf = true;
                itemSize = RPTREE_NODE_LEAF_ITEM_SIZE;

                RPTreeLeafNode* tempP= new RPTreeLeafNode();
                thisNode =  tempP ;
            } else {
                isLeaf = false;
                itemSize = RPTREE_NODE_CHILD_ITEM_SIZE;
                thisNode = new RPTreeChildNode();
            }
            //nodeCount++;

            uint16_t itemCount;
            byte empty_b;



            fis->read( reinterpret_cast<char*>(&empty_b) , sizeof(byte) );
            empty_b=endian::LittleByte(empty_b);

            fis->read( reinterpret_cast<char*>(&itemCount) , sizeof(uint16_t) );
            itemCount=endian::LittleShort(itemCount);

          //  uint32_t itemBlockSize = itemCount * itemSize;

          //  buffer = new byte[itemBlockSize];            // allocate buffer for item sisze

            // get the node items - leaves or child nodes
            uint32_t startChromID, endChromID;
            uint32_t startBase, endBase;
            int itemBlockSize = itemCount * itemSize;

            //Readitems into buffer and set to stringstream
            std::vector<char>  buffer(itemBlockSize);
            fis->read(&buffer[0],itemBlockSize);
            std::stringstream       itemStream;
            //Set string stream, does not copy buffer
            itemStream.rdbuf()->pubsetbuf(&buffer[0],itemBlockSize);

            //TODO make item reading work with buffer
            for (uint32_t item = 0; item < itemCount; ++item) {

                // always extract the bounding rectangle

                //Read from buffer
                itemStream.read( reinterpret_cast<char*>(&startChromID) , sizeof(uint32_t) );
                startChromID=endian::LittleLong(startChromID);

                itemStream.read( reinterpret_cast<char*>(&startBase) , sizeof(uint32_t) );
                startBase=endian::LittleLong(startBase);

                itemStream.read( reinterpret_cast<char*>(&endChromID) , sizeof(uint32_t) );
                endChromID=endian::LittleLong(endChromID);

                itemStream.read( reinterpret_cast<char*>(&endBase) , sizeof(uint32_t) );
                endBase=endian::LittleLong(endBase);

                if (isLeaf) {
                    uint64_t dataOffset;
                    uint64_t dataSize;


                    itemStream.read( reinterpret_cast<char*>(&dataOffset) , sizeof(uint64_t) );
                    dataOffset=endian::LittleDouble(dataOffset);

                    itemStream.read( reinterpret_cast<char*>(&dataSize) , sizeof(uint64_t) );
                    dataSize=endian::LittleDouble(dataSize);
                    //TODO Check creation here, it was incorrect.
                    thisNode->insertItem(new RPTreeLeafNodeItem(startChromID, startBase, endChromID, endBase,
                            dataOffset, dataSize));
                } else {
                    // get the child node pointed to in the node item

                    uint64_t nodeOffset=0;

                    itemStream.read( reinterpret_cast<char*>(&nodeOffset) , sizeof(uint64_t) );
                    nodeOffset=endian::LittleDouble(nodeOffset);


                    // Recursive call to read next child node
                    // The test on chromIds is designed to stop the descent when the tree reaches the level of an
                    // individual chromosome.  These are loaded later on demand.
                    RPTreeNode* childNode;
                    if (startChromID != endChromID || forceDescend) {
                        childNode = readRPTreeNode(fis, nodeOffset, forceDescend);
                    } else {
                        childNode = new RPTreeNodeProxy(*fis, nodeOffset, startChromID);
                    }

                    // insert child node item
                    thisNode->insertItem(new RPTreeChildNodeItem(startChromID, startBase, endChromID,
                            endBase, childNode));

                }
                fileOffset += itemSize;
            }

        } catch (...) {
            throw new std::runtime_error("Error reading R+ tree nodes: \n");
        }

        // return success
        return thisNode;
    }
