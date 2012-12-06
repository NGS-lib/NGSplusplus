#ifndef UWRITERGFF_H_INCLUDED
#define UWRITERGFF_H_INCLUDED



#include <iostream>
#include <fstream>
#include "uWriterBase.h"
namespace NGS {

class uWriterGFF : public uWriterBase {
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterGFF() {}
	void writeToken(const uToken& token);
private:
	static DerivedRegister<uWriterGFF> reg;
    bool hasStarted=false;

}; // End of class uWriterBed4

} // End of namespace NGS
#endif // UWRITERGFF_H_INCLUDED
