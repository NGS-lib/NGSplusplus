#ifndef UWRITERGTF_H_INCLUDED
#define UWRITERGTF_H_INCLUDED

#include <iostream>
#include <fstream>
#include "uWriterBase.h"
namespace NGS {

class uWriterGTF : public uWriterBase {
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterGTF() {}
	void writeToken(const uToken& token);
private:
	static DerivedRegister<uWriterGTF> reg;
    bool hasStarted=false;

};
} // End of namespace NGS


#endif // UWRITERGFF_H_INCLUDED

