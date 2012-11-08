#ifndef UWRITEWIG_H_INCLUDED
#define UWRITEWIG_H_INCLUDED
#include <iostream>
#include <fstream>
#include "uWriterBase.h"
namespace NGS {

class uWriterWig : public uWriterBase {
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterWig() {}
	void writeToken(const uToken& token);

private:
	static DerivedRegister<uWriterWig> reg;
    bool chrValidate(const std::string& chr);
    int span=0;


}; // End of class uWriterBed4

} // End of namespace NGS
#endif // UWRITEWIG_H_INCLUDED
