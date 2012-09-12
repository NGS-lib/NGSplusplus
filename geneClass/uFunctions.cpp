#include "uTags.h"
#include "uRegion.h"
#include "uFunctions.h"
namespace NGS {
//Some functions using Regions, we can put elsewhere later.

namespace RegionTags{

    int getTagCount(std::string chr, int start, int stop, uTagsExperiment* ourExp, OverlapType overlap)
    {
    if (ourExp->isModeGradual())
        {
            std::cerr << "ERROR in function getTagDensity: File open in progressive read" <<std::endl;
            abort();
        }
       return (ourExp->getSubsetCount(chr,start, stop,overlap));
    }

    int getTagCount(uRegion inputRegion, uTagsExperiment* ourExp, OverlapType overlap)
    {
       return getTagCount(inputRegion.getChr(), inputRegion.getStart(), inputRegion.getEnd(),ourExp,overlap ); //     .setDensity(ourExp->getRegionCount(returnRegion.getChr(),returnRegion.getStart(), returnRegion.getEnd(), false));

    }

    void getSignal(uRegion* inputRegion, uTagsExperiment* ourExp, bool overlap){

    if (ourExp->isModeGradual())
        {
            std::cerr << "ERROR in function getTagDensity: File open in progressive read" <<std::endl;
            abort();
        }

       inputRegion->setSignal(ourExp->getRegionSignal(inputRegion->getChr(),inputRegion->getStart(), inputRegion->getEnd(), overlap));

    }


}
} // End of namespace NGS
