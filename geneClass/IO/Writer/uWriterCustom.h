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
	uWriterCustom();
	virtual void writeToken(const uToken& token);
	void setFieldsNames(const std::vector<std::string>& fieldsNames);
	void setDelimiter(char delimiter) { m_delimiter = delimiter; }

private:
	static DerivedRegister<uWriterCustom> reg;

	std::vector<std::string> m_fieldsNames;
	char m_delimiter = '\t';

}; // End of class uWriterCustom

} // End of namespace NGS
#endif // UWRITERCUSTOM_H_INCLUDED
