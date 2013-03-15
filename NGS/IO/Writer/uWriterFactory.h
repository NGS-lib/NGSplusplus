
#ifndef UWRITERFACTORY_H_INCLUDED
#define UWRITERFACTORY_H_INCLUDED

#include "uWriterBase.h"
#include "uWriterBed4.h"
#include "uWriterBed6.h"
#include "uWriterBedGraph.h"
#include "uWriterCustom.h"
#include "uWriterGFF.h"
#include "uWriterGTF.h"
#include "uWriterSam.h"
#include "uWriterWig.h"
namespace NGS{

template<typename T> uWriterBase * createT() { return new T; }

struct uWriterBaseFactory {
    typedef std::map<std::string, std::function<uWriterBase*()> > Writer_map_type;
    virtual ~uWriterBaseFactory(){};
    static std::shared_ptr<uWriterBase> createInstance(std::string const& s) {
        Writer_map_type::iterator it = getWriterMap()->find(s);

        if(it == getWriterMap()->end()){
            throw uWriter_exception_base()<<string_error("Asked for unregistered type: "+s+" in Writer, failling");
        }
      //  std::cerr << "Returning" <<std::endl;
        return std::shared_ptr<uWriterBase>(it->second());
    }
    static uWriterBaseFactory * GetFact()
        {
        static uWriterBaseFactory instance;
        return &instance;
        }

    void Register(const std::string &s, std::function<uWriterBase*()> pCreatFunc)
    {
            auto it = getWriterMap()->find(s);
            if(it == getWriterMap()->end())
            {
                 getWriterMap()->insert(std::pair<std::string, std::function<uWriterBase*() >> (s, pCreatFunc));
            }
            else
            {
                throw uWriter_exception_base()<<string_error("Duplicated type registering in Writer, failling\n Type is:"+s+"\n");
            }
    }


    uWriterBaseFactory()
    {
        Register("BED", &uWriterBed4::Create);
        Register("BED4", &uWriterBed4::Create);
        Register("BED6", &uWriterBed6::Create);
        Register("WIG", &uWriterWig::Create);
        Register("SAM", &uWriterSam::Create);
        Register("GFF", &uWriterGFF::Create);
        Register("GTF", &uWriterGTF::Create);
        Register("BEDGRAPH", &uWriterBedGraph::Create);
        Register("CUSTOM", &uWriterCustom::Create);

    }
protected:
    static Writer_map_type * mapItem;
    static Writer_map_type * getWriterMap()
    {
        if(!mapItem)
	{
	    mapItem= new Writer_map_type;
	}
        return mapItem;
    };

private:


};
//
//template<typename T>
//struct DerivedWriterRegister : uWriterBaseFactory
//{
//    ~DerivedWriterRegister(){};
//    DerivedWriterRegister(std::string const& s)
//    {
//        auto it = getWriterMap()->find(s);
//        if(it == getWriterMap()->end())
//	{
//             getWriterMap()->insert(std::pair<std::string, std::function<uWriterBase*() >> (s, &createT<T>));
//        }
//	else
//        {
//            throw uWriter_exception_base()<<string_error("Duplicated type registering in Writer, failling\n Type is:"+s+"\n");
//        }
//    }
//};

}
#endif
