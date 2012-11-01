#include "uParserBase.h"

namespace NGS 
{

/** \brief Destructor.
 */
uParserBase::~uParserBase() 
{
    if (m_dynamicStream == true) 
    {
        delete m_pIostream;
    }
    m_pIostream = NULL;
}

/** \brief Initialize the uParserBase object (open file and parse header).
 * \param const std::string& filename: name of the bed file to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserBase::init(const std::string& filename, bool header) 
{
    std::ifstream* ifs = new std::ifstream(filename.c_str(), std::ifstream::in);
    if (!ifs->is_open())
    {
        std::string error = "Error opening file: " + filename;
        throw std::runtime_error(error.c_str());
    }
    else
    {
        m_pIostream = ifs;
        m_dynamicStream = true;
    }
}

/** \brief Initialize the uParserBase object (set stream and parse header).
 * \param std::iostream* stream: name of the bed stream to parse.
 * \param bool header: true if there is a header to parse (value at false by default).
 */
void uParserBase::init(std::iostream* stream, bool header) 
{
    m_pIostream = stream;
    m_dynamicStream = false;
}


/** \brief Default constructor.
 */
uParserBase::uParserBase(){};

/** \brief Check if we are at the end of the file/stream
 * \return true is we are at the end, otherwise false.
 */
bool uParserBase::eof() const
{ 
    return m_pIostream->peek() == EOF; 
}

std::map<std::string, std::function<uParserBase*()> > *uParserBaseFactory::mapItem;
}
