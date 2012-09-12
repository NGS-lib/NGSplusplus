#ifndef IPARSERBASE_H_INCLUDED
#define IPARSERBASE_H_INCLUDED

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

class uToken;
class header_param;
class parserBase{
    public :
    parserBase();
    virtual ~parserBase();

     virtual void init(const std::string& filename, bool header = false)=0;
	 virtual void init(std::iostream* stream, bool header = false)=0;
	 virtual void init(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t')=0;
	 virtual void init(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t')=0;
	parserBase& operator=(const parserBase& copyFrom) = delete;
	parserBase(const parserBase&) = delete;
	/** \brief Check if input data is at end of file.
	  */
	bool eof() const { return m_pIostream->peek() == EOF; }
	virtual uToken getNextEntry();
	/** \brief Get a specific data from header.
	  */
	//std::string getHeaderParam(header_param name) const { return m_headerData.getParam(name); }
	/** \brief Check if there is a value associated with a given param.
	  * \param header_param& name: name of the param to check.
	  */
//	bool isHeaderParamSet(const header_param& name) const { return m_headerData.isParamSet(name); }
private:
    std::istream* m_pIostream=nullptr;
};

//Thank you Stack Overflow for this basic structure

template<typename T> parserBase * createT() { return new T; }

struct parserBaseFactory {
   typedef std::map<std::string, std::function<parserBase*()> > map_type;


    static std::shared_ptr<parserBase> createInstance(std::string const& s) {
        map_type::iterator it = getMap()->find(s);
        if(it == getMap()->end())
            return nullptr;
        return std::shared_ptr<parserBase>(it->second());
    }

protected:
 static map_type * mapItem;
        static map_type * getMap() {
        if(!mapItem) { mapItem= new map_type; }
        return mapItem;
    };

};
template<typename T>
struct DerivedRegister : parserBaseFactory {
    DerivedRegister(std::string const& s) {
       // getMap()->insert(std::pair(s, &createT<T>));
       map_type * test= getMap();
      (*test)[s]= std::bind(&createT<T>);
    }
};
#endif // IPARSERBASE_H_INCLUDED
