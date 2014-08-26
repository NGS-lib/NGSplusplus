
// ***************************************************************************
//   BBFileReader.cpp (c) 2014
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



//TODO Refactor everything to stop using raw pointers. C++ 11 here we come!


#include "BBFileReader.h"
#include <vector>
#include <string>

    BBFileReader::BBFileReader(){
    }

    BBFileReader::~BBFileReader(){

           if (fileHeader_!=NULL) {
                delete(fileHeader_);
           }

           if (zoomLevels_!=NULL) {
                delete(zoomLevels_);
           }

           if (totalSummaryBlock_!=NULL) {
                delete(totalSummaryBlock_);
           }

           if (chromosomeIDTree_!=NULL) {
                delete(chromosomeIDTree_);
           }

           if (chromosomeDataTree_!=NULL) {
                delete(chromosomeDataTree_);
           }

    }


    BBFileReader::BBFileReader(std::string path, std::ifstream& stream) {
        open(path,stream);
    }

    void BBFileReader::open(std::string path, std::ifstream& stream){

        //Get endian information
        endian::InitEndian();
        fis_ = &stream;

        // read in file header
        fileOffset_ = BBFILE_HEADER_OFFSET;
        if (fis_->good()==false)
            throw std::runtime_error("Invalid stream passed to BBFileReader ");

        fileHeader_ = new BBFileHeader(path, *fis_, fileOffset_);
        //fileHeader_.print();

        if (!fileHeader_->isHeaderOK()) {
            throw std::runtime_error("Error reading BBFile header for: ");
        }

        // get data characteristics
        uncompressBufSize_ = fileHeader_->getUncompressBuffSize();

        // update file offset past BBFile header
        fileOffset_ += BBFileHeader::BBFILE_HEADER_SIZE;

        // get zoom level count from file header
        zoomLevelCount_ = fileHeader_->getZoomLevels();

        // extract zoom level headers and zoom data records
        // Note: zoom headers Table D immediately follow the BBFile Header
        if (zoomLevelCount_ > 0) {

            zoomLevelOffset_ = fileOffset_;

            zoomLevels_ = new BBZoomLevels(*fis_, zoomLevelOffset_, zoomLevelCount_,
                     uncompressBufSize_);

            // end of zoom level headers - compare with next BBFile item location
            fileOffset_ += zoomLevelCount_ * BBZoomLevelHeader::ZOOM_LEVEL_HEADER_SIZE;
        }

        // get the AutoSQL custom BigBed fields
        autoSqlOffset_ = fileHeader_->getAutoSqlOffset();
        if (autoSqlOffset_ != 0) {
            // read in .as entry
            // mFileOffset = mAutoSqlOffset + sizeof(.as format field);
        }

        // get the Total Summary Block (Table DD)
        fileOffset_ = fileHeader_->getTotalSummaryOffset();
        if (fileHeader_->getVersion() >= 2 && fileOffset_ > 0) {
            totalSummaryBlock_ = new BBTotalSummaryBlock(*fis_, fileOffset_);
            fileOffset_ += BBTotalSummaryBlock::TOTAL_SUMMARY_BLOCK_SIZE;
        }

        // get Chromosome Data B+ Tree (Table E, F, G, H) : should always exist
        chromIDTreeOffset_ = fileHeader_->getChromosomeTreeOffset();
        if (chromIDTreeOffset_ != 0) {
            fileOffset_ = chromIDTreeOffset_;
            chromosomeIDTree_ = new BPTree(*fis_, fileOffset_);
        }

        // get R+ chromosome data location tree (Tables K, L, M, N)
        chromDataTreeOffset_ = fileHeader_->getFullIndexOffset();
        if (chromDataTreeOffset_ != 0) {
            fileOffset_ = chromDataTreeOffset_;
            bool forceDescend = false;
            chromosomeDataTree_ = new RPTree(*fis_, fileOffset_, uncompressBufSize_, forceDescend);
        }

        // get number of data records indexed by the R+ chromosome data location tree
        fileOffset_ = fileHeader_->getFullDataOffset();
        dataCount_ = getDataCount(*fis_, fileOffset_);
}



    /*
    *   Method returns the Big Binary File pathname.
    *
    *   Returns:
    *       Big Binary File pathname
    * */

    std::string BBFileReader::getBBFilePath() {
        return path_;
    }

    /*
    *   Method returns the Big Binary File input stream handle.
    *
    *   Returns:
    *       Big Binary File input stream handle
    * */

    std::ifstream* BBFileReader::getBBFis() {
        return fis_;
    }

    /*
    *   Method returns the Big Binary File header which identifies
    *   the file type and content.
    *
    *   Returns:
    *       Big Binary File header (Table C)
    * */

    BBFileHeader* BBFileReader::getBBFileHeader() {
        return fileHeader_;
    }

    /*
    *   Method returns if the Big Binary File is BigBed.
    *
    *   Returns:
    *       Boolean identifies if Big Binary File is BigBed
    *       (recognized from magic number in file header Table C)
    * */

    bool BBFileReader::isBigBedFile() {
        return fileHeader_->isBigBed();
    }

    /*
    *   Method returns if the Big Binary File is BigWig
    *
    *   Returns:
    *       Boolean identifies if Big Binary File is BigWig
    *       (recognized from magic number in file header Table C)
    * */

    bool BBFileReader::isBigWigFile() {
        return fileHeader_->isBigWig();
    }

    /*
    *   Method returns the total number of data records in the file.
    *
    *   Returns:
    *       Count of the total number of compressed/uncompressed data records:
    *           which for BigBed is the number of bed features,
    *           and for BiGWifg is the number of wig sections.
    * */

    uint32_t BBFileReader::getDataCount() {
        return dataCount_;
    }

    /*
    *   Method returns the number of chromosomes/contigs in the file.
    *
    *   Note: This is itemCount from B+ tree header BBFile Table E.
    *
    *   Returns:
    *       Count of the total number of chromosomes/contigs in the file.
    * */

    int64_t BBFileReader::getChromosomeNameCount() {
        return chromosomeIDTree_->getItemCount();
    }

    /*
    *   Method returns the number of chromosome/contig regions in the file.
    *
    *   Note: This is itemCount from R+ tree header BBFile Table K.
    *
    *   Returns:
    *       Count of the total number of chromosome/contig regions in the file.
    * */

    int64_t BBFileReader::getChromosomeRegionCount() {
        return chromosomeDataTree_->getItemCount();
    }

    /*
   *   Method returns the Big Binary File decompressed buffer size.
   *
   *   Returns:
   *       Largest required buffer size for decompressed data chunks (from Table C)
   * */

    int32_t BBFileReader::getDecompressionBufSize() {
        return uncompressBufSize_;
    }


    /*
    *   Method returns the total summary block for the Big Binary File.
    *
    *   Returns:
    *       Total summary block data statistics for the whole file (Table DD)
    * */

    BBTotalSummaryBlock* BBFileReader::getTotalSummaryBlock() {
        return totalSummaryBlock_;
    }

    /*
    *   Method returns the B+ Chromosome Index Tree.
    *
    *   Returns:
    *       B+ Chromosome Index Tree (includes Tables  E, F, G, H)
    * */

    BPTree* BBFileReader::getChromosomeIDTree() {
        return chromosomeIDTree_;
    }

    /*
    *   Method returns the R+ Chromosome Data Locations Tree.
    *
    *   Returns:
    *       R+ Chromosome Data Locations Tree (includes Tables  K, L, M, N)
    * */

    RPTree* BBFileReader::getChromosomeDataTree() {
        return chromosomeDataTree_;
    }

    /*
    *   Method returns number of zoom level data is included in the file.
    *
    *   Returns:
    *       Number of zoom levels (from Table C)
    * */

    int32_t BBFileReader::getZoomLevelCount() {
        return zoomLevelCount_;
    }

    /*
    *   Method returns the zoom levels in the Big Binary File.
    *
    *   Returns:
    *       Zoom level object containing zoom level headers and R+ zoom data locations tree
    *       (includes Tables  D, O)
    * */

    BBZoomLevels* BBFileReader::getZoomLevels() {
        return zoomLevels_;
    }

    /*
    *   Method finds the zoom data bounds in R+ tree for a chromosome ID range.
    *
    *   Parameters:
    *       zoomLevel - zoom level
    *       startChromID - start chromosome for the region
    *       endChromID - end chromosome for the region
    *
    *   Returns:
    *       Chromosome region bounds for chromosome ID range
    * */

    RPChromosomeRegion* BBFileReader::getZoomLevelBounds(int32_t zoomLevel, int32_t startChromID,
                                                 int32_t endChromID) {

        RPChromosomeRegion* chromosomeBounds =  zoomLevels_->getZoomLevelRPTree(zoomLevel)->getChromosomeRegion(startChromID, endChromID);

        return chromosomeBounds;
    }

    /*
    *   Method finds chromosome bounds for entire chromosome ID range in the zoom level R+ tree.
    *
    *   Parameters:
    *       zoomLevel - zoom level
    *
    *   Returns:
    *       Chromosome bounds for the entire chromosome ID range in the R+ tree.
    * */

    RPChromosomeRegion* BBFileReader::getZoomLevelBounds(int32_t zoomLevel) {

        RPChromosomeRegion* chromosomeBounds =
                zoomLevels_->getZoomLevelRPTree(zoomLevel)->getChromosomeBounds();

        return chromosomeBounds;
    }

    /*
    *   Method returns the zoom record count for the zoom level.
    *
    *   Parameters:
    *       zoomLevel - zoom level
    *
    *   Returns:
    *       Chromosome bounds for the entire chromosome ID range in the R+ tree.
    * */

    int32_t BBFileReader::getZoomLevelRecordCount(int32_t zoomLevel) {

        return zoomLevels_->getZoomLevelFormats().at(zoomLevel - 1)->getZoomRecordCount();
    }

    /*
    *   Method finds chromosome key name for the associated chromosome ID in the B+ tree.
    *
    *   Returns:
    *       chromosome key name for associated chromosome ID.
    * */

     std::string BBFileReader::getChromosomeName(int32_t chromID) {

        std::string chromosomeName = chromosomeIDTree_->getChromosomeName(chromID);
        return chromosomeName;
    }

    /*
    *   Method finds chromosome names in the B+ chromosome index tree.
    *
    *   Returns:
    *       LIst of all chromosome key names in the B+ tree.
    * */

    std::vector<std::string> BBFileReader::getChromosomeNames() {

        std::vector<std::string> chromosomeList = chromosomeIDTree_->getChromosomeNames();
        return chromosomeList;
    }

    /*
   *   Returns a chromosome ID  which  can be used to search for a
   *   corresponding data section in the R+ tree for data.
   *
      Parameters:
   *       chromosomeKey - chromosome name of valid key size.
   *
   *
   *   Note: A chromosomeID of -1 means chromosome name not included in B+ tree.
   *
   * */

    int32_t BBFileReader::getChromosomeID(std::string chromosomeKey) {

        int32_t chromosomeID = chromosomeIDTree_->getChromosomeID(chromosomeKey);

        return chromosomeID;
    }

    /*
    *   Method finds the chromosome bounding region in R+ tree for a chromosome ID range.
    *
    *   Parameters:
    *       startChromID - start chromosome for the region
    *       endChromID - end chromosome for the region
    *
    *   Returns:
    *       Chromosome region bounds for chromosome ID range
    * */

    RPChromosomeRegion* BBFileReader::getChromosomeBounds(int32_t startChromID, int32_t endChromID) {

        RPChromosomeRegion* chromosomeBounds =
                chromosomeDataTree_->getChromosomeRegion(startChromID, endChromID);

        return chromosomeBounds;
    }

    /*
    *   Method finds the chromosome bounds for the entire chromosome ID range in the R+ tree.
    *
    *   Returns:
    *       chromosome bounds for the entire chromosome ID range in the R+ tree.
    * */

    RPChromosomeRegion* BBFileReader::getChromosomeBounds() {

        RPChromosomeRegion* chromosomeBounds = chromosomeDataTree_->getChromosomeBounds();

        return chromosomeBounds;
    }

    /*
    *   Method finds all chromosome data regions in the R+ tree.
    *
    *   Returns:
    *       List of chromosome ID's and regions.
    * */

     std::vector<RPChromosomeRegion*> BBFileReader::getChromosomeRegions() {

        std::vector<RPChromosomeRegion*> regionList = chromosomeDataTree_->getAllChromosomeRegions();

        return regionList;
    }

    /*
    *   Method finds all zoom level data regions in the R+ tree.
    *
    *   Parameters:
    *       int zoomLevel - zoom level
    *   Returns:
    *       List of chromosome ID's and regions for the zoom level.
    * */

    std::vector<RPChromosomeRegion*> BBFileReader::getZoomLevelRegions(int32_t zoomLevel) {

        std::vector<RPChromosomeRegion*> regionList =
                zoomLevels_->getZoomLevelRPTree(zoomLevel)->getAllChromosomeRegions();

        return regionList;
    }

//    /**
//     * Returns an iterator for BigBed features which occupy a chromosome selection region.
//     * <p/>
//     * Note: the BBFile type should be BigBed; else a null iterator is returned.
//     * <p/>
//     * Parameters:
//     * startChromosome - name of start chromosome
//     * startBase     - starting base position for features
//     * endChromosome - name of end chromosome
//     * endBase       - ending base position for feature
//     * contained     - flag specifies bed features must be contained in the specified
//     * base region if true; else can intersect the region if false
//     * <p/>
//     * Returns:
//     * Iterator to provide BedFeature(s) for the requested chromosome region.
//     * Error conditions:
//     * 1) An empty iterator is returned if region has no data available
//     * 2) A null object is returned if the file is not BigBed.(see isBigBedFile method)
//     */
//    BigBedIterator* BBFileReader::getBigBedIterator(std:: startChromosome, int32_t startBase,
//                                            std:: endChromosome, int32_t endBase, bool contained) {
//
//        if (!isBigBedFile())
//            return NULL;
//
//
//        // go from chromosome names to chromosome ID region
//        RPChromosomeRegion* selectionRegion = getChromosomeBounds(startChromosome, startBase,
//                endChromosome, endBase);
//
//        // check for valid selection region
//        if (selectionRegion == NULL)
//            return new BigBedIterator();  // an empty iterator
//
//        // compose an iterator
//        BigBedIterator* bedIterator = new BigBedIterator(fis_, chromosomeIDTree_, chromosomeDataTree_,
//                selectionRegion_, contained);
//
//        return bedIterator;
//    }

//    /**
//     * Returns an iterator for BigBed features for all chromosome regions.
//     * <p/>
//     * Note: the BBFile type should be BigBed; else a null iterator is returned.
//     * <p/>
//     * Returns:
//     * Iterator to provide BedFeature(s) for all chromosome regions.
//     * Error conditions:
//     * 1) An empty iterator is returned if region has no data available
//     * 2) A null object is returned if the file is not BigBed.(see isBigBedFile method)
//     */
//    BigBedIterator* BBFileReader::getBigBedIterator() {
//
//        if (!isBigBedFile())
//            return NULL;
//
//
//        // get all region bounds
//        RPChromosomeRegion* selectionRegion = chromosomeDataTree_->getChromosomeBounds();
//
//        // compose an iterator
//        bool contained = true;   /// all regions are contained
//        BigBedIterator* bedIterator = new BigBedIterator(fis_, chromosomeIDTree_, chromosomeDataTree_,
//                selectionRegion_, contained);
//
//        return bedIterator;
//    }


    /**
     * Returns an iterator for BigWig values which occupy the specified startChromosome region.
     * <p/>
     * Note: the BBFile type should be BigWig; else a null iterator is returned.
     * <p/>
     * Parameters:
     * startChromosome  - name of start chromosome
     * startBase    - starting base position for features
     * endChromosome  - name of end chromosome
     * endBase      - ending base position for feature
     * contained    - flag specifies bed features must be contained in the specified
     * base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     * 2) A null object is returned if the file is not BigWig.(see isBigWigFile method)
     */
     BigWigIterator BBFileReader::getBigWigIterator(std::string startChromosome, int32_t startBase,
                                            std::string endChromosome, int32_t endBase, bool contained) {


        if (!isBigWigFile())
            throw std::runtime_error("Cannot provide iterator");

        // go from chromosome names to chromosome ID region
        RPChromosomeRegion* selectionRegion = getChromosomeBounds(startChromosome, startBase,
                endChromosome, endBase);

        // check for valid selection region, return empty iterator if null
        if (selectionRegion == NULL)
            return BigWigIterator();

        // compose an iterator
      //  BigWigIterator* wigIterator = new ;

        return BigWigIterator(fis_, chromosomeIDTree_, chromosomeDataTree_,
                                selectionRegion, contained);
    }

    /**
     * Returns an iterator for BigWig values for all chromosome regions.
     * <p/>
     * Note: the BBFile type should be BigWig; else a null iterator is returned.
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for all chromosome regions.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     * 2) A null object is returned if the file is not BigWig.(see isBigWigFile method)
     */
    BigWigIterator BBFileReader::getBigWigIterator() {
    try {
        if (!isBigWigFile())
            throw std::runtime_error("Could not provide iterator");

        // get all regions bounds
        RPChromosomeRegion* selectionRegion = chromosomeDataTree_->getChromosomeBounds();

        // compose an iterator
        bool contained = true;       // all regions are containe

        return BigWigIterator(fis_, chromosomeIDTree_, chromosomeDataTree_,
                selectionRegion, contained);

     }
     catch(...){
     throw;}
    }


    /**
     * Returns an iterator for zoom level records for the chromosome selection region.
     * <p/>
     * Note: the BBFile can be BigBed or BigWig.
     * <p/>
     * Parameters:
     * zoomLevel - zoom level for data extraction; levels start at 1
     * startChromosome - start chromosome name
     * startBase     - staring base position for features
     * endChromosome - end chromosome name
     * endBase       - ending base position for feature
     * contained     - flag specifies bed features must be contained in the
     * specified base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     */
    ZoomLevelIterator BBFileReader::getZoomLevelIterator(unsigned int zoomLevel, std::string startChromosome, int startBase,
                                                  std::string endChromosome, int endBase, bool contained) {
        // check for valid zoom level
        if (zoomLevel < 1 || (unsigned int)zoomLevel > zoomLevelCount_)
            throw std::runtime_error("Error: ZoomLevelIterator zoom level is out of range\n");

        // get the appropriate zoom level R+ zoom data index tree
        RPTree* zoomDataTree = zoomLevels_->getZoomLevelRPTree(zoomLevel);

        // go from chromosome names to chromosome ID region
        RPChromosomeRegion* selectionRegion = getChromosomeBounds(startChromosome, startBase,
                endChromosome, endBase);

        // check for valid selection region
        if (selectionRegion == NULL) {
            return EmptyIterator();
        }

        return ZoomLevelIterator(fis_, chromosomeIDTree_,
                zoomDataTree, zoomLevel, selectionRegion, contained);
    }

    /**
     * Returns an iterator for zoom level records for all chromosome regions.
     * <p/>
     * Note: the BBFile can be BigBed or BigWig.
     * <p/>
     * Parameters:
     * zoomLevel - zoom level for data extraction; levels start at 1
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     */
    ZoomLevelIterator BBFileReader::getZoomLevelIterator(unsigned int zoomLevel) {

        // check for valid zoom level
        if (zoomLevel < 1 || zoomLevel > zoomLevelCount_)
            throw std::runtime_error("Error: ZoomLevelIterator zoom level is out of range\n");

        // get the appropriate zoom level R+ zoom data index tree
        RPTree* zoomDataTree = zoomLevels_->getZoomLevelRPTree(zoomLevel);

        // get all regions bounds
        RPChromosomeRegion* selectionRegion = zoomDataTree->getChromosomeBounds();

        // compose an iterator
        bool contained = true;   //all regions are contained

        return ZoomLevelIterator(fis_, chromosomeIDTree_,
                zoomDataTree, zoomLevel, selectionRegion, contained);
    }

    /**
     * Returns an iterator for zoom level records for the chromosome selection region.
     * <p/>
     * Note: the BBFile can be BigBed or BigWig.
     * <p/>
     * Parameters:
     * zoomLevel - zoom level for data extraction; levels start at 1
     * selectionRegion - chromosome selection region consists of:
     * startChromID - ID of starting chromosome
     * startBase     - staring base position for features
     * endChromID - ID of endind chromosome
     * endBase       - ending base position for feature
     * contained     - flag specifies bed features must be contained in the
     * specified base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     */
    ZoomLevelIterator BBFileReader::getZoomLevelIterator(unsigned int zoomLevel, RPChromosomeRegion* selectionRegion,
                                                  bool contained) {
        // check for valid zoom level
        if (zoomLevel < 1 || (unsigned int)zoomLevel > zoomLevelCount_)
            throw new std::runtime_error("Error: ZoomLevelIterator zoom level is out of range\n");

        // get the appropriate zoom level R+ zoom data index tree
        RPTree* zoomDataTree = zoomLevels_->getZoomLevelRPTree(zoomLevel);

        return ZoomLevelIterator(fis_, chromosomeIDTree_,zoomDataTree, (int32_t)zoomLevel, selectionRegion, contained);
    }

    /*
    *   Method generates a chromosome bounds region for the supplied chromosome region name.
    *
    *   Note: No attempt is made to verify the region exists in the file data, nor
    *   which data is being examined.
    *
    *   Parameters:
    *       startChromosome - name of start chromosome
    *       startBase - starting base position for region
    *       endChromosome - name of end chromosome
    *       endBase - ending base position for region
    *
    *   Returns:
    *       Chromosome bounds of a named chromosome region for data extraction;
    *       or null for regions not found in the B+ chromosome index tree.
    * */

    RPChromosomeRegion* BBFileReader::getChromosomeBounds(std::string startChromosome, int startBase,
                                                   std::string endChromosome, int endBase) {

        // If the chromosome name length is > the key size we can't distinguish it
        if (startChromosome.length() > chromosomeIDTree_->getKeySize()) {
            return NULL;
        }

        // find the chromosome ID's using the name to get a valid name key, then associated ID
        std::string startChromKey = chromosomeIDTree_->getChromosomeKey(startChromosome);
        int startChromID = chromosomeIDTree_->getChromosomeID(startChromKey);
        if (startChromID < 0)       // mChromosome not in data?
            return NULL;

        std::string endChromKey = chromosomeIDTree_->getChromosomeKey(endChromosome);
        int endChromID = chromosomeIDTree_->getChromosomeID(endChromKey);
        if (endChromID < 0)       // mChromosome not in data?
            return NULL;

        // create the bounding mChromosome region
        RPChromosomeRegion* chromBounds = new RPChromosomeRegion(startChromID, startBase,
                endChromID, endBase);

        return chromBounds;
    }

    /*
    *   Method reads data count which heads the data section of the BBFile.
    *
    *   Returns:
    *       Data count of the number of data records:
    *          number of Bed features for BigBed
    *          number of Wig sections for BigWig
    * */

    uint32_t BBFileReader::getDataCount(std::ifstream& fis_, uint64_t fileOffset) {
        int dataCount;
   //     LittleEndianInputStream lbdis = null;
  //      DataInputStream bdis = null;

        // Note: dataCount in BBFile is simply a 4 byte int
        // positioned at fullDataOffset in Table C
    //    char[] buffer = new byte[4];

        try {
            // read dataCount into a buffer
            fis_.clear();
            fis_.seekg(fileOffset);
             if (fis_.eof())
                std::cerr<<"Hit end of file in seekg in getDataCount\n";
      //      fis_.readFully(buffer);
             //   lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
                fis_.read( reinterpret_cast<char*>(&dataCount) , sizeof(int32_t) );
                dataCount=endian::LittleLong(dataCount);

        } catch (...) {
            throw new std::runtime_error("Error reading data count for all data");
        }

        // data count was read properly
        return dataCount;
    }
