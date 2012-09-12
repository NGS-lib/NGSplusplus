#ifndef IPARSERBASE_H_INCLUDED
#define IPARSERBASE_H_INCLUDED


class parserBase{
    public :
    parserBase();
    virtual ~parserBase();

     virtual void init(const std::string& filename, bool header = false);
	 virtual void init(std::iostream* stream, bool header = false);
	 virtual void init(const std::string& filename, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	 virtual void init(std::iostream* stream, const std::vector<std::string>& fieldsNames, bool header = false, char delimiter = '\t');
	~parserBase();
	parserBase& operator=(const parserBase& copyFrom) = delete;
	parserBase(const parserBase&) = delete;
	/** \brief Check if input data is at end of file.
	  */
	bool eof() const { return m_pIostream->peek() == EOF; }
	virtual uToken getNextEntry();
	/** \brief Get a specific data from header.
	  */
	std::string getHeaderParam(header_param name) const { return m_headerData.getParam(name); }
	/** \brief Check if there is a value associated with a given param.
	  * \param header_param& name: name of the param to check.
	  */
	bool isHeaderParamSet(const header_param& name) const { return m_headerData.isParamSet(name); }
private:
    istream* m_pIostream;
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
#endif // IPARSER_H_INCLUDED
#endif // IPARSERBASE_H_INCLUDED
