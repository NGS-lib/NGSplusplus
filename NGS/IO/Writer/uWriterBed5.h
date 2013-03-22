#ifndef UWRITERBED5_H_INCLUDED
#define UWRITERBED5_H_INCLUDED

#include "uWriterBed.h"
#include <iostream>
#include <fstream>

namespace NGS
{

class uWriterBed5 : public uWriterBed
{
public:
    /** \brief Empty constructor (call object with init through factory instead)
      */
    uWriterBed5() {}
    virtual void writeToken(const uToken& token);
    static uWriterBase * Create() { return new uWriterBed5(); }
private:

}; // End of class uWriterBed5

} // End of namespace NGS
#endif // UWRITERBED5_H_INCLUDED
