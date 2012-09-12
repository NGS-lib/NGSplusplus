#ifndef UWRITERBED_H_INCLUDED
#define UWRITERBED_H_INCLUDED

#include <iostream>
#include "../NGS++.h"
#include "uWriterBase.h"

namespace NGS {

class uWriterBed : public uWriterBase {
public:
	uWriterBed() {}
	~uWriterBed() {} 
	void init(const std::string& filename);
	void init(std::ofstream* os);
	void writeToken(const uToken& token);

private:
	ostream* m_pOstream;
	static DerivedRegister<uWriterBed> reg;

}; // End of class uWriterBed

} // End of namespace NGS
#endif // UWRITERBED_H_INCLUDED
