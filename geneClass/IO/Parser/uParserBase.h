#ifndef IuParserBase_H_INCLUDED
#define IuParserBase_H_INCLUDED

#include <map>
#include <string>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "../uToken.h"
#include "../../uGeneException.h"
#include "../uHeader.h"
namespace NGS {

class uParserBase{
    public :
    uParserBase();
    virtual ~uParserBase(){};

     virtual void init(const std::string& filename, bool header = false)=0;
	 virtual void init(std::iostream* stream, bool header = false)=0;
	 virtual void init(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t')=0;
	 virtual void init(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t')=0;
	uParserBase& operator=(const uParserBase& copyFrom) = delete;
	uParserBase(const uParserBase&) = delete;
	/** \brief Check if input data is at end of file.
	  */
	bool eof() const;
	virtual uToken getNextEntry()=0;
	/** \brief Get an unformated version of header (i.e.: a single string containing the whole header)
	*/
	std::string getUnformatedHeader() const { return m_headerData.getUnformatedHeader(); }
	/** \brief Get a specific data from header.
	  */
	std::string getHeaderParam(header_param name) const { return m_headerData.getParam(name); } ;
    std::vector<std::string> getHeaderParamVector(header_param name) const{return m_headerData.getParamVector(name); } ;
	/** \brief Check if there is a value associated with a given param.
	  * \param header_param& name: name of the param to check.
	  */
	bool isHeaderParamSet(const header_param& name) const { return m_headerData.isParamSet(name); }

protected:
	bool m_dynamicStream = false;
    std::istream* m_pIostream=nullptr;
    uHeader m_headerData;

};

//Thank you Stack Overflow for this basic structure

template<typename T> uParserBase * createT() { return new T; }

struct uParserBaseFactory {
   typedef std::map<std::string, std::function<uParserBase*()> > parser_map_type;
    virtual ~uParserBaseFactory(){};
    static std::shared_ptr<uParserBase> createInstance(std::string const& s) {
        parser_map_type::iterator it = getParserMap()->find(s);

      //  std::cout << getParserMap()->size() <<std::endl;
       //  parser_map_type::iterator temp= getParserMap()->begin();
      //   std::cout << it->first <<std::endl;
        if(it == getParserMap()->end()){
            throw uParser_exception_base()<<string_error("Asked for unregistered type, failling");
        }
      //  std::cerr << "Returning" <<std::endl;
        return std::shared_ptr<uParserBase>(it->second());
    }

protected:
 static parser_map_type * mapItem;
        static parser_map_type * getParserMap() {
        if(!mapItem) { mapItem= new parser_map_type; }
        return mapItem;
    };

};
template<typename T>
struct DerivedParserRegister : uParserBaseFactory {
    ~DerivedParserRegister(){};
    DerivedParserRegister(std::string const& s) {

         auto it = getParserMap()->find(s);
        if(it == getParserMap()->end()){
             getParserMap()->insert(std::pair<std::string, std::function<uParserBase*() >> (s, &createT<T>));
        } else
        {
            throw uParser_exception_base()<<string_error("Duplicated type registering in Parser, failling\n Type is:"+s+"\n");
        }


    }

};
}
#endif // IuParserBase_H_INCLUDED
