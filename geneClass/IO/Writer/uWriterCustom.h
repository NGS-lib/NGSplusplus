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
//	void init(const std::string& filename, const std::vector<std::string>& fieldsNames);
//	void init(std::ostream* os, const std::vector<std::string>& fieldsNames);
//	void setFieldsNames(const std::vector<std::string> fieldsNames);
	virtual void writeToken(const uToken& token);

private:
//	std::vector<std::string> m_fieldsNames;
	static DerivedRegister<uWriterCustom> reg;

}; // End of class uWriterCustom

} // End of namespace NGS
#endif // UWRITERCUSTOM_H_INCLUDED
