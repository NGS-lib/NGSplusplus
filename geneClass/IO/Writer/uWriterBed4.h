#ifndef UWRITERBED4_H_INCLUDED
#define UWRITERBED4_H_INCLUDED

#include <iostream>
#include <fstream>
#include "../NGS++.h"
#include "uWriterBed.h"

namespace NGS {

class uWriterBed4 : public uWriterBed {
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterBed4() {}
	virtual ~uWriterBed4() {}
	virtual void init(const std::string& filename) { uWriterBed::init(filename); }
	virtual void init(std::ostream* os) { uWriterBed::init(os); }
	virtual void writeToken(const uToken& token);

private:
	static DerivedRegister<uWriterBed4> reg;

}; // End of class uWriterBed4

} // End of namespace NGS
#endif // UWRITERBED4_H_INCLUDED
