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
    void setValueType(const token_param type);
private:
	static DerivedRegister<uWriterWig> reg;
    bool chrValidate(const std::string& chr);
    int m_span=0;
    std::string m_chr="";
    token_param m_scoreType=token_param::SCORE;

}; // End of class uWriterBed4

} // End of namespace NGS
#endif // UWRITEWIG_H_INCLUDED
