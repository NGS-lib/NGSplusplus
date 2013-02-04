#ifndef UWRITERBED4_H_INCLUDED
#define UWRITERBED4_H_INCLUDED

#include "uWriterBed.h"
#include <iostream>
#include <fstream>

namespace NGS
{

class uWriterBed4 : public uWriterBed
{
public:
    /** \brief Empty constructor (call object with init through factory instead)
      */
    uWriterBed4() {}
    virtual void writeToken(const uToken& token);
     static uWriterBase * Create() { return new uWriterBed4(); }
private:
   // static DerivedRegister<uWriterBed4> reg;

}; // End of class uWriterBed4

} // End of namespace NGS
#endif // UWRITERBED4_H_INCLUDED
