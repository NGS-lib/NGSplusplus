#include "WigItem.h"
#include "RPChromosomeRegion.h"
#include "RPTreeLeafNodeItem.h"
#include "BigWigSectionHeader.h"

#include <stdlib.h>
#include <stdexcept>
#include <map>
#include <fstream>
#include <sstream>
#include "endian_helper.h"
using namespace std;


//This should replace the use of a BigWigSection, as it seems to be created then deleted in the BigWigDataBlock section

    int32_t extractSectionData(RPChromosomeRegion& selectionRegion, std::vector<char>& sectionBuffer,
                           std::map<uint32_t, std::string>& chromosomeMap,bool contained,
                            std::vector<WigItem>& wigItemList) {

        std::stringstream myStream;
        myStream.rdbuf()->pubsetbuf(&sectionBuffer[0],sectionBuffer.size());
        BigWigSectionHeader wigSectionHeader(myStream);


        //wigSectionHeader_ = new BigWigSectionHeader(myStream_);


        // check for valid Wig item type
        if(wigSectionHeader.getItemType() == WigTypeNamespace::Unknown)
        {
            throw std::runtime_error("Read error on wig section leaf index ");
        }

        // include header in data segment size accounting
        int sectionDataSize = wigSectionHeader.SECTION_HEADER_SIZE;


        // get the section's data item specifications
        // Note: A RuntimeException is thrown if wig section header is not read properly
        int32_t chromID =  wigSectionHeader.getChromID();
        std::string chromosome = (chromosomeMap)[chromID];
        int32_t itemCount = wigSectionHeader.getItemCount();
        int32_t chromStart = wigSectionHeader.getChromosomeStart();
        int32_t chromEnd = wigSectionHeader.getChromosomeEnd();
        int32_t itemStep = wigSectionHeader.getItemStep();
        int32_t itemSpan =  wigSectionHeader.getItemSpan();
        int32_t itemIndex = 0;
        int32_t startBase = 0;
        int32_t endBase = 0;
        float value = 0.0f;

        // find Wig data type - BBFile Table J item type
        WigTypeNamespace::WigItemType itemType = wigSectionHeader.getItemType();

        // check if all leaf items are selection hits
        RPChromosomeRegion* itemRegion = new RPChromosomeRegion(chromID, chromStart,
                            chromID, chromEnd);
        int32_t leafHitValue = itemRegion->compareRegions(&selectionRegion);

        // extract Wig data records
        // Note: the buffer input stream is positioned past section header
        try {
            for(int32_t index = 0; index < itemCount; ++index) {
                ++itemIndex;

                    if(itemType == WigTypeNamespace::FixedStep){
                        startBase = chromStart;
                        endBase = startBase + itemSpan;

                        myStream.read( reinterpret_cast<char*>(&value) , sizeof(float) );
                        value=endian::LittleFloat(value);
                        chromStart = startBase + itemStep;
                        sectionDataSize += BigWigSectionHeader::FIXEDSTEP_ITEM_SIZE;
                    }
                    else if(itemType == WigTypeNamespace::VarStep){


                        myStream.read( reinterpret_cast<char*>(&startBase) , sizeof(uint32_t) );
                        startBase=endian::LittleLong(startBase);

                        endBase = startBase + itemSpan;

                        myStream.read( reinterpret_cast<char*>(&value) , sizeof(float) );
                        value=endian::LittleFloat(value);
                        sectionDataSize += BigWigSectionHeader::VARSTEP_ITEM_SIZE;

                    }
                    else if(itemType == WigTypeNamespace::BedGraph){

                        myStream.read( reinterpret_cast<char*>(&startBase) , sizeof(uint32_t) );
                        startBase=endian::LittleLong(startBase);

                        myStream.read( reinterpret_cast<char*>(&endBase) , sizeof(uint32_t) );
                        endBase=endian::LittleLong(endBase);

                        myStream.read( reinterpret_cast<char*>(&value) , sizeof(float) );
                        value=endian::LittleFloat(value);

                        sectionDataSize += BigWigSectionHeader::BEDGRAPH_ITEM_SIZE;
                    }

                // contained leaf region items are always added - otherwise test conditions
                if(leafHitValue == 0) {
                   // WigItem* bbItem = new ;
                    wigItemList.push_back(WigItem(itemIndex, chromosome, startBase, endBase, value));
                }
                else {
                    itemRegion = new RPChromosomeRegion(chromID, startBase, chromID, endBase);
                    int32_t itemHitValue = itemRegion->compareRegions(&selectionRegion);


                    //TODO: Stop working with pointers? Yes?
                    // hitValue < 2 needed for intersection; hitValue < 1 needed for contained = true
                    if(itemHitValue == 0 || ( !contained && abs(itemHitValue) < 2) ) {


                       // WigItem* bbItem = new ;
                        wigItemList.push_back(WigItem(itemIndex, chromosome, startBase, endBase, value));
                    }
                }

            }

        }catch(...) {
            throw new std::runtime_error("Read error for Wig section item ");
        }

        return sectionDataSize;
    }
