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
	~uWriterBed();
	void init(const std::string& filename);
	void init(std::ostream* os);
	virtual void writeToken(const uToken& token);

}; // End of class uWriterBed

} // End of namespace NGS
#endif // UWRITERBED_H_INCLUDED
