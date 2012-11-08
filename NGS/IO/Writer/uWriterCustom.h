#ifndef UWRITERCUSTOM_H_INCLUDED
#define UWRITERCUSTOM_H_INCLUDED

#include "uWriterBase.h"
#include <iostream>
#include <fstream>
namespace NGS
{

class uWriterCustom : public uWriterBase
{
public:
    /** \brief Empty constructor (call object with init through factory instead)
      */
    uWriterCustom();
    virtual ~uWriterCustom() {};
    virtual void writeToken(const uToken& token);

private:
    static DerivedRegister<uWriterCustom> reg;

}; // End of class uWriterCustom

} // End of namespace NGS
#endif // UWRITERCUSTOM_H_INCLUDED
