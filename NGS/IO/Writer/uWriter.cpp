#include "uWriter.h"
#include "uWriterBase.h"
#include "uWriterFactory.h"
namespace NGS
{

/** \brief Constructor with filename.
  * \param const std::string& filename: name of the file to write output.
  * \param const std::string& type: Type of file to write output (i.e: Bed, Sam, etc...).
  */
uWriter::uWriter(const std::string& filename, const std::string& type)
{
   // uWriterBaseFactory myFactory;
   // m_pWriterBase = myFactory.createInstance(type);
    try
    {
        m_pWriterBase=uWriterBaseFactory::GetFact()->createInstance(type);
        m_pWriterBase->init(filename);
    }
    catch (std::runtime_error& e)
    {
        throw e;
    }
}

/** \brief Constructor with stream.
  * \param std::ostream* os: stream to write output.
  * \param const std::string& type: Type of file to write output (i.e: Bed, Sam, etc...).
  */
uWriter::uWriter(std::ostream* os, const std::string& type)
{
  //  uWriterBaseFactory myFactory;
   // m_pWriterBase = myFactory.createInstance(type);
    m_pWriterBase=uWriterBaseFactory::GetFact()->createInstance(type);
    m_pWriterBase->init(os);
}

/** \brief Custom constructor with filename.
  * \param const std::string& filename: name of the file to write output.
  * \param const std::vector<std::string>& fieldsNames: Vector containing the name of every fields (columns). In correct order.
  * \param const std::string& type: Type of file to write output (i.e: Bed, Sam, etc...).
  */
uWriter::uWriter(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter)
{
   // uWriterBaseFactory myFactory;
   // m_pWriterBase = myFactory.createInstance("CUSTOM");
    try
    {
        m_pWriterBase=uWriterBaseFactory::GetFact()->createInstance("CUSTOM");
        m_pWriterBase->init(filename);
    }
    catch (std::runtime_error& e)
    {
        throw e;
    }
    m_pWriterBase->setFieldsNames(fieldsNames);
    m_pWriterBase->setDelimiter(delimiter);
}

/** \brief Custom constructor with stream.
  * \param std::ostream* os: stream to write output.
  * \param const std::vector<std::string>& fieldsNames: Vector containing the name of every fields (columns). In correct order.
  * \param const std::string& type: Type of file to write output (i.e: Bed, Sam, etc...).
  */
uWriter::uWriter(std::ostream* os, const std::vector<std::string>& fieldsNames, char delimiter)
{
   // uWriterBaseFactory myFactory;
   // m_pWriterBase = myFactory.createInstance("CUSTOM");
    m_pWriterBase=uWriterBaseFactory::GetFact()->createInstance("CUSTOM");
    m_pWriterBase->init(os);
    m_pWriterBase->setFieldsNames(fieldsNames);
    m_pWriterBase->setDelimiter(delimiter);
}

/** \brief Print the information of the token in the desired format.
  * \param const uToken& token: The token to print.
  */
void uWriter::writeToken(const uToken& token)
{
    m_pWriterBase->writeToken(token);
}

/** \brief Print an unformated string.
  * \param const std::string& str: the string to write to the file or stream.
  */
void uWriter::printString(const std::string& str)
{
    m_pWriterBase->printString(str);
}

/** \brief setHeaderParameter to be writen
 *
 * \param param header_param
 * \param value std::string
 * \return void
 *
 */
void uWriter::addToHeader(header_param param,std::string value){
    m_pWriterBase->addToHeader(param, value);
}

/** \brief Write the header, ideally use this once
 *
 * \param param header_param
 * \param value std::string
 * \return void
 *
 */
void uWriter::writeHeader(){
    m_pWriterBase->writeHeader();
}

} // End of namespace NGS
