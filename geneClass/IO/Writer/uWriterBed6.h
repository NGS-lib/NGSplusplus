#ifndef UWRITERBED6_H_INCLUDED
#define UWRITERBED6_H_INCLUDED

#include <iostream>
#include <fstream>
#include "../NGS++.h"
#include "uWriterBed.h"

namespace NGS {

class uWriterBed6 : public uWriterBed {
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterBed6() {}
	virtual ~uWriterBed6() {}
	virtual void init(const std::string& filename) { uWriterBed::init(filename); }
	virtual void init(std::ostream* os) { uWriterBed::init(os); }
	virtual void writeToken(const uToken& token);

private:
	static DerivedRegister<uWriterBed6> reg;

}; // End of class uWriterBed6

} // End of namespace NGS
#endif // UWRITERBED6_H_INCLUDED
