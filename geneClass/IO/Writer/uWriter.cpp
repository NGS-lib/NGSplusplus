#include "uWriter.h"

namespace NGS {

/** \brief Constructor with filename. 
  * \param const std::string& filename: name of the file to write output.
  * \param const std::string& type: Type of file to write output (i.e: Bed, Sam, etc...). 
  */
uWriter::uWriter(const std::string& filename, const std::string& type) {
	uWriterBaseFactory myFactory;
	m_pWriterBase = myFactory.createInstance(type);
	m_pWriterBase->init(filename);
}

/** \brief Constructor with stream.
  * \param std::ostream* os: stream to write output.
  * \param const std::string& type: Type of file to write output (i.e: Bed, Sam, etc...). 
  */
uWriter::uWriter(std::ostream* os, const std::string& type) {
	uWriterBaseFactory myFactory;
	m_pWriterBase = myFactory.createInstance(type);
	m_pWriterBase->init(os);
}

/** \brief Print the information of the token in the desired format.
  * \param const uToken& token: The token to print.
  */
void uWriter::writeToken(const uToken& token) {
	m_pWriterBase->writeToken(token);
}

}; // End of namespace NGS
