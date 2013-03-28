// ***************************************************************************
// uWriterBed.h (c) 2013
// Alexei Nordell-Markovits : Sherbrooke University
// Charles Joly Beauparlant : Laval University
//
//       This file is part of the NGS++ library.
//
//    The NGS++ library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program (lgpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// ***************************************************************************

#ifndef UWRITERBED_H_INCLUDED
#define UWRITERBED_H_INCLUDED

#include "uWriterBase.h"
#include <iostream>
#include <fstream>

namespace NGS
{

class uWriterBed : public uWriterBase
{
public:
    /** \brief Empty constructor (call object with init through factory instead)
      */
    uWriterBed() {}
    /** \brief Destructor.
     */
    virtual ~uWriterBed() {}
    virtual void writeToken(const uToken& token);

}; // End of class uWriterBed

} // End of namespace NGS
#endif // UWRITERBED_H_INCLUDED
