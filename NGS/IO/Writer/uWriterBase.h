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
#include "../../utility/utility.h"
namespace NGS
{


namespace UCSCHeader
{
const char UCSCCOMMENT='#';
const std::string UCSCBROWSER="browser";
const std::string UCSCPOS="position";
const std::string UCSCHIDE="hide";
const std::string UCSCHIDEALL="hide all";
const std::string UCSCDENSE="dense";
const std::string UCSCDENSEALL="dense all";
const std::string UCSCPACK="pack";
const std::string UCSCPACKALL="pack all";
const std::string UCSCSQUISH="squish";
const std::string UCSCSQUISHALL="squish all";
const std::string UCSCFULLALL="full all";
const std::string UCSCFULL="full";


/**< All definition */

template <typename ...Tail>
std::string buildBrowserLine(const std::string& pType, const std::string& curTrack, Tail&&... tail)
{
    std::string concatValues= utility::STRING::concatStringListWithSpaces(curTrack,std::forward<Tail>(tail)...);
    return (UCSCBROWSER+" "+pType+" "+concatValues+"\n");
}

enum class UCSCBrowseType
{
    HIDE,PACK,SQUISH,FULL,DENSE
};
enum class UCSCBrowseTypeMult
{
    HIDE_ALL,DENSE_ALL,PACK_ALL,SQUISH_ALL,FULL_ALL
};

/**< Return line to write to specify "position" tag for UCSC browser */
static std::string getUCSCPositionLine(const std::string & pChr, long long pStart, long long pEnd)
{
    return (UCSCBROWSER+" "+UCSCPOS+" "+pChr+":"+std::to_string(pStart)+"-"+std::to_string(pEnd)+"\n");
}

static std::string getUCSCLine(const UCSCBrowseType& pType)
{
    switch (pType)
    {
    case UCSCBrowseType::PACK:
        return UCSCBROWSER+" "+UCSCPACK+"\n";
    case UCSCBrowseType::HIDE:
        return UCSCBROWSER+" "+UCSCHIDE+"\n";
    case UCSCBrowseType::SQUISH:
        return UCSCBROWSER+" "+UCSCSQUISH+"\n";
    case UCSCBrowseType::FULL:
        return UCSCBROWSER+" "+UCSCFULL+"\n";
    case UCSCBrowseType::DENSE:
        return UCSCBROWSER+" "+UCSCDENSE+"\n";
    }
     /**< Cannot reach here */
    assert(false);
}

template <typename ...Tail>
std::string getUCSCLine(const UCSCBrowseTypeMult& pType,const std::string& curTrack, Tail&&... tail)
{
    /**< Valid with list */
    switch (pType)
    {
    case UCSCBrowseTypeMult::PACK_ALL:
        return buildBrowserLine(UCSCPACKALL,curTrack,tail...);
    case UCSCBrowseTypeMult::HIDE_ALL:
        return buildBrowserLine(UCSCHIDEALL,curTrack,tail...);
    case UCSCBrowseTypeMult::SQUISH_ALL:
        return buildBrowserLine(UCSCSQUISHALL,curTrack,tail...);
    case UCSCBrowseTypeMult::FULL_ALL:
        return buildBrowserLine(UCSCFULLALL,curTrack,tail...);
    case UCSCBrowseTypeMult::DENSE_ALL:
        return buildBrowserLine(UCSCDENSEALL,curTrack,tail...);
    }
    /**< Cannot reach here */
   assert(false);
}

}
class uWriterBase
{
public :
    uWriterBase() {}
    virtual ~uWriterBase();
    virtual void init(const std::string& filename);
    virtual void init(std::ostream* os);
    virtual void writeToken(const uToken& token)=0;
    void printString(const std::string& str);
    uWriterBase& operator=(const uWriterBase& copyFrom) = delete;
    uWriterBase(const uWriterBase&) = delete;
    void setFieldsNames(const std::vector<std::string>& fieldsNames);
    void setDelimiter(char delimiter)
    {
        m_delimiter = delimiter;
    }

    virtual void addToHeader(header_param param,std::string value);
    virtual void writeHeader();

protected:
    std::ostream* m_pOstream = nullptr;
    bool m_dynamicStream = false;
    std::vector<std::string> m_fieldsNames= {};
    char m_delimiter = '\t';
    uHeader m_headerData=uHeader();

};

//Thank you Stack Overflow for this basic structure

//template<typename T> uWriterBase * createT()
//{
//    return new T;
//}
//
//struct uWriterBaseFactory
//{
//    typedef std::map<std::string, std::function<uWriterBase*()> > map_type;
//    virtual ~uWriterBaseFactory() {}
//    static std::shared_ptr<uWriterBase> createInstance(std::string const& s)
//    {
//        map_type::iterator it = getMap()->find(s);
//        if(it == getMap()->end())
//    {
//
//        throw uWriter_invalid_type_instance()<<string_error("Asked for unregistered type: "+s+" in Writer, failling");
//
//        }
//        return std::shared_ptr<uWriterBase>(it->second());
//    }
//
//protected:
//    static map_type* mapItem;
//        static map_type* getMap()
//    {
//        if(!mapItem)
//    {
//            mapItem= new map_type;
//        }
//        return mapItem;
//    };
//}; // End of struct uWriterBaseFactory
//
//template<typename T>
//struct DerivedRegister : uWriterBaseFactory
//{
//    ~DerivedRegister(){};
//    DerivedRegister(std::string const& s)
//    {
//        getMap()->insert(std::pair<std::string, std::function<uWriterBase*() >> (s, &createT<T>));
//    }
//}; // End of struct DerivedRegister

} // End of namespace NGS
#endif // UWRITERBASE_H_INCLUDED
