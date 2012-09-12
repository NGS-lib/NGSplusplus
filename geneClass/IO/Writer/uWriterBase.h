#ifndef UWRITERBASE_H_INCLUDED
#define UWRITERBASE_H_INCLUDED
#include <iostream>

namespace NGS {

class uWriterBase{
public :
	uWriterBase() {}
	virtual ~uWriterBase() {}
	void init(const std::string& filename)=0;
	void init(std::ostream* os)=0;
	void writeToken(const uToken& token);
	uWriterBase& operator=(const uParserBase& copyFrom) = delete;
	uWriterBase(const uParserBase&) = delete;
private:
	std::ostream* m_pOstream = nullptr;
};

//Thank you Stack Overflow for this basic structure

template<typename T> uWriterBase * createT() { return new T; }

struct uWriterBaseFactory {
	typedef std::map<std::string, std::function<uWriterBase*()> > map_type;

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
	DerivedRegister(std::string const& s) {
		map_type * test= getMap();
		(*test)[s]= std::bind(&createT<T>);
	}
}; // End of struct DerivedRegister 

} // End of namespace NGS
#endif // UWRITERBASE_H_INCLUDED
