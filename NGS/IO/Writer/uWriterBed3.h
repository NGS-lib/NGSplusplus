#ifndef UWRITERBED3_H_INCLUDED
#define UWRITERBED3_H_INCLUDED

#include "uWriterBed.h"
#include <iostream>
#include <fstream>

namespace NGS
{

class uWriterBed3 : public uWriterBed
{
public:
    /** \brief Empty constructor (call object with init through factory instead)
      */
    uWriterBed3() {}
    virtual void writeToken(const uToken& token);
    static uWriterBase * Create() { return new uWriterBed3(); }
private:

}; // End of class uWriterBed3

} // End of namespace NGS
#endif // UWRITERBED3_H_INCLUDED
