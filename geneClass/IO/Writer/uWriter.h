#ifndef UWRITER_H_INCLUDED
#define UWRITER_H_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "../../NGS++.h"
//#include "../../uGeneException.h"
#include "uWriterBase.h"
//#include "../uToken.h"
//#include "../uHeader.h"

namespace NGS {

class header;

class uWriter {
public:
	uWriter(const std::string& filename, const std::string& type);
	uWriter(std::ostream* os, const std::string& type);
	uWriter(const std::string& filename, const std::vector<std::string>& fieldsNames, const std::string& type);
	uWriter(std::ostream* os, const std::vector<std::string>& fieldsNames, const std::string& type);
	void writeToken(const uToken& token);

private:
	uHeader m_headerData;
	std::shared_ptr<uWriterBase> m_pWriterBase = nullptr;
}; // End of class uWriter

} // End of namespace NGS
#endif // UWRITER_H_INCLUDED
