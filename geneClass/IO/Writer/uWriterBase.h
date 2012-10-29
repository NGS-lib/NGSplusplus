#ifndef UWRITERBASE_H_INCLUDED
#define UWRITERBASE_H_INCLUDED
#include <iostream>
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

class uWriterBase{
public :
	uWriterBase() {}
	virtual ~uWriterBase();
	void init(const std::string& filename);
	void init(std::ostream* os);
	virtual void writeToken(const uToken& token)=0;
	void printString(const std::string& str);
	uWriterBase& operator=(const uWriterBase& copyFrom) = delete;
	uWriterBase(const uWriterBase&) = delete;
	void setFieldsNames(const std::vector<std::string>& fieldsNames);
	void setDelimiter(char delimiter) { m_delimiter = delimiter; }

	virtual void addToHeader(header_param param,std::string value);
	virtual void writeHeader();

protected:
	std::ostream* m_pOstream = nullptr;
	bool m_dynamicStream = false;
	std::vector<std::string> m_fieldsNames={};
	char m_delimiter = '\t';
	uHeader m_headerData;
};

//Thank you Stack Overflow for this basic structure

template<typename T> uWriterBase * createT() { return new T; }

struct uWriterBaseFactory {
	typedef std::map<std::string, std::function<uWriterBase*()> > map_type;
	virtual ~uWriterBaseFactory() {}
	static std::shared_ptr<uWriterBase> createInstance(std::string const& s) {
		map_type::iterator it = getMap()->find(s);
		if(it == getMap()->end()) {
			return nullptr;
		}
		return std::shared_ptr<uWriterBase>(it->second());
	}

protected:
	static map_type* mapItem;
        static map_type* getMap() {
		if(!mapItem) {
			mapItem= new map_type;
		}
		return mapItem;
	};
}; // End of struct uWriterBaseFactory

template<typename T>
struct DerivedRegister : uWriterBaseFactory {
	~DerivedRegister(){};
	DerivedRegister(std::string const& s) {
		getMap()->insert(std::pair<std::string, std::function<uWriterBase*() >> (s, &createT<T>));
	}
}; // End of struct DerivedRegister

} // End of namespace NGS
#endif // UWRITERBASE_H_INCLUDED
