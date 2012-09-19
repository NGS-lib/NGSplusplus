#ifndef UWRITERCUSTOM_H_INCLUDED
#define UWRITERCUSTOM_H_INCLUDED

#include <iostream>
#include <fstream>
#include "../NGS++.h"
#include "uWriterBase.h"

namespace NGS {

class uWriterCustom : public uWriterBase { 
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterCustom() {}
	~uWriterCustom();
	virtual void init(const std::string& filename);
	virtual void init(std::ostream* os);
	virtual void writeToken(const uToken& token);

private:
	static DerivedRegister<uWriterCustom> reg;

}; // End of class uWriterCustom

} // End of namespace NGS
#endif // UWRITERCUSTOM_H_INCLUDED
