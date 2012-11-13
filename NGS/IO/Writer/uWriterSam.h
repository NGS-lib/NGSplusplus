#ifndef UWRITERSAM_H_INCLUDED
#define UWRITERSAM_H_INCLUDED

#include <iostream>
#include <fstream>
#include "uWriterBase.h"
namespace NGS {

class uWriterSam : public uWriterBase {
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterSam() {}
	void writeToken(const uToken& token);
    void writeHeader();
private:
	static DerivedRegister<uWriterSam> reg;
    bool chrValidate(const std::string& chr);
    bool hasStarted=false;


}; // End of class uWriterBed4

} // End of namespace NGS
#endif

