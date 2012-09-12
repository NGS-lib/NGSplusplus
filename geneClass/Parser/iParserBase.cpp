#include "base.h"
std::map<std::string, std::function<parserBase*()> > *parserBaseFactory::mapItem;

parserBase::parserBase(){};

parserBase::~parserBase(){};
