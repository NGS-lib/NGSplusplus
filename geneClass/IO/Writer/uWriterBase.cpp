#include "uWriterBase.h"

namespace NGS
{

/** \brief Destructor
  */
uWriterBase::~uWriterBase()
{
    if (m_dynamicStream == true)
    {
        delete m_pOstream;
    }
    m_pOstream = nullptr;
}

/** \brief Initialise the Bed writer with a file name.
  * \param const std::string& filename: name of the file to write output.
  */
void uWriterBase::init(const std::string& filename) {
	if (filename.size() == 0) {
		throw std::runtime_error("Filename must be longer than 0");
	}
	std::filebuf fb;
	fb.open (filename.c_str(),std::ios::out);
	if (!fb.is_open()) {
		std::string error = "Error opening file: " + filename;
		throw std::runtime_error(error.c_str());
	}
	else {
		std::ostream* os = new std::ostream(&fb);
		m_pOstream = os;
	}
	m_dynamicStream = true;
}

/** \brief Initialise the Bed writer with a stream
  * \param std::ostream* os: the stream to save the data to.
  */
void uWriterBase::init(std::ostream* os)
{
    if (os != nullptr && os->good() == true)
    {
        m_pOstream = os;
    }
    else
    {
        throw std::runtime_error("Invalid stream.");
    }
    m_dynamicStream = false;
}

/** \brief Print an unformated string.
  * \param const std::string& str: the string to write to the file or stream.
  */
void uWriterBase::printString(const std::string& str)
{
    *m_pOstream << str;
}

/** \brief Set the fields name (only used for the Custom file)
  */
void uWriterBase::setFieldsNames(const std::vector<std::string>& fieldsNames)
{
    if (fieldsNames.size() == 0)
    {
        throw no_fields_names() << string_error("fieldsNames vector is empty");
    }
    m_fieldsNames = fieldsNames;
}


void uWriterBase::addToHeader(header_param param,std::string value){
    m_headerData._addToParam(param, value);
}

void uWriterBase::writeHeader(){

}

std::map<std::string, std::function<uWriterBase*()> > *uWriterBaseFactory::mapItem;

} // End of namespace NGS
