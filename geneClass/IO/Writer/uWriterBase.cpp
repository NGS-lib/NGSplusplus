#include "uWriterBase.h"

namespace NGS {

void uWriterBase::setFieldsNames(const std::vector<std::string>& fieldsNames) {
	if (fieldsNames.size() == 0) {
		throw no_fields_names() << string_error("fieldsNames vector is empty");
	}
	m_fieldsNames = fieldsNames;
}

void uWriterBase::printString(const std::string& str) {
	*m_pOstream << str;
}

std::map<std::string, std::function<uWriterBase*()> > *uWriterBaseFactory::mapItem;

} // End of namespace NGS
