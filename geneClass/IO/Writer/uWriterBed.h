#ifndef UWRITERBED_H_INCLUDED
#define UWRITERBED_H_INCLUDED

#include <iostream>
#include <fstream>
#include "../NGS++.h"
#include "uWriterBase.h"

namespace NGS {

class uWriterBed : public uWriterBase { 
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterBed() {}
	virtual void writeToken(const uToken& token);

}; // End of class uWriterBed

} // End of namespace NGS
#endif // UWRITERBED_H_INCLUDED
