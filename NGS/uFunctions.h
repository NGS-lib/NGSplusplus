#ifndef UFUNCTIONS_H_INCLUDED
#define UFUNCTIONS_H_INCLUDED
namespace NGS {

namespace RegionTags{
//How many tags exactly fit in the region
int getTagCount(std::string chr, int start, int stop, uTagsExperiment* ourExp, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
int getTagCount(uRegion inputRegion, uTagsExperiment* ourExp, OverlapType overlap=OverlapType::OVERLAP_PARTIAL);
void getSignal(uRegion* inputRegion, uTagsExperiment* ourExp,  bool signal);
}

} // End of namespace NGS

#endif

