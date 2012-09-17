#include "uWriterBase.h"

namespace NGS {

std::map<std::string, std::function<uWriterBase*()> > *uWriterBaseFactory::mapItem;

} // End of namespace NGS
