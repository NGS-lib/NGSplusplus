#ifndef BBZOOMLEVELS_H
#define BBZOOMLEVELS_H


// ***************************************************************************
//   BBZoomLevels.h (c) 2014
//   Copyright @ Alexei Nordell-Markovits : Sherbrooke University
//
//    This file is part BWReader library.
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


#include <string>
#include <vector>
#include "stdint.h"
#include "BBZoomLevelHeader.h"
#include "BBZoomLevelFormat.h"
#include "RPTree.h"

class BBZoomLevels
{
    public:
        BBZoomLevels();
        virtual ~BBZoomLevels();

         BBZoomLevels(std::ifstream& fis, uint64_t fileOffset, uint32_t zoomLevels
                        , uint32_t uncompressBufSize);
        uint64_t getZoomHeadersOffset();
        uint32_t getZoomHeaderCount();
        std::vector<BBZoomLevelHeader*> getZoomLevelHeaders();
        BBZoomLevelHeader* getZoomLevelHeader(uint32_t level);
        std::vector<BBZoomLevelFormat*> getZoomLevelFormats();
        RPTree* getZoomLevelRPTree(uint32_t level);
        uint32_t readZoomHeaders(std::ifstream& fis, uint64_t fileOffset, uint32_t zoomLevels);
    protected:
    private:
        uint64_t zoomHeadersOffset_;
        uint32_t zoomLevelsCount_;
        std::vector<BBZoomLevelHeader*> zoomLevelHeaders_;
        std::vector<BBZoomLevelFormat*> zoomLevelFormatList_;
        std::vector<RPTree*> zoomLevelRPTree_;
};

#endif // BBZOOMLEVELS_H
