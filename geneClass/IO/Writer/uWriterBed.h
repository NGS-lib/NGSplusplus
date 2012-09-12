#ifndef UWRITERBED_H_INCLUDED
#define UWRITERBED_H_INCLUDED

#include <iostream>
#include "../NGS++.h"

namespace NGS {

class uWriterBed : public uWriterBase {
public:
	uWriterBed(const std::string& filename);
	uWriterBed(std::ostream* os);
	void writeToken(const uToken& token);

private:
	ostream* m_pOstream;

}; // End of class uWriterBed

} // End of namespace NGS
#endif // UWRITERBED_H_INCLUDED
