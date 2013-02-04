#include "../../uGeneException.h"
#include "uParser.h"
#include "uParserBase.h"
//#include "uParserBed.h"
#include "uParserFactory.h"
namespace NGS
{

/** \brief Filename default constructor.
 * \param const std::string& filename: the name of the file to parse.
 * \param const std::string& type: the type of file (i.e.: "BED")
 * \param bool header: true if there is a header (value at false by default).
 */
uParser::uParser(const std::string& filename, const std::string & type, bool header)
{

   // uParserBaseFactory myFact;
  //  DerivedParserRegister<uParserBed> uParserBed::reg("BED");

    try
    {
         m_pParserBase=uParserBaseFactory::GetFact()->createInstance(type);
     //   m_pParserBase=myFact.createInstance(type);
        m_pParserBase->init(filename, header);
    }
    catch(...)
    {
        throw;
        //  std::cerr <<fetchStringError(e);
    }
}

/** \brief Stream default constructor.
 * \param std::iostream* stream: the stream to parse.
 * \param const std::string& type: the type of stream (i.e.: "BED")
 * \param bool header: true if there is a header (value at false by default).
 */
uParser::uParser(std::istream* stream, const std::string & type, bool header)
{

  //  uParserBaseFactory myFact;
  m_pParserBase= uParserBaseFactory::GetFact()->createInstance(type);
   // m_pParserBase=myFact.createInstance(type);
    m_pParserBase->init(stream, header);

}

/** \brief Filename custom constructor.
 * \param const std::string& filename: the name of the file to parse.
 * \param const std::vector<std::string>& fieldsNames: The name of every column in the file.
 * \param bool header: true if there is a header (value at false by default).
 */
uParser::uParser(const std::string& filename, const std::vector<std::string>& fieldsNames, char delimiter)
{

   // uParserBaseFactory myFact;
    m_pParserBase=uParserBaseFactory::GetFact()->createInstance("CUSTOM");
  //  m_pParserBase=myFact.createInstance("CUSTOM");
    m_pParserBase->init(filename, fieldsNames, delimiter);

}

/** \brief Stream custom constructor.
 * \param std::iostream* stream: The stream to parse.
 * \param const std::vector<std::string>& fieldsNames: The name of every column in the file.
 * \param bool header: true if there is a header (value at false by default).
 */
uParser::uParser(std::istream* stream, const std::vector<std::string>& fieldsNames, char delimiter)
{

    //uParserBaseFactory myFact;
     m_pParserBase=uParserBaseFactory::GetFact()->createInstance("CUSTOM");
   // m_pParserBase=myFact.createInstance("CUSTOM");
    m_pParserBase->init(stream, fieldsNames, delimiter);

}

/** \brief Destructor.
 */
uParser::~uParser() {}

/** \brief Check if we are at the end of the file (or of the stream)
 * \return true if we are at the end of the file, otherwise return false.
 */
bool uParser::eof() const
{
    return m_pParserBase->eof();
}

/** \brief Create a token from current point in the file (or the stream).
 * \return a uToken objet containing the infos for the next entry.
 */
uToken uParser::getNextEntry()
{
    return  m_pParserBase->getNextEntry();
}

}
