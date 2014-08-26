
/**< Unfortunately endian.h is not standard */



/*
Copyright (C) 2003  Promit Roy

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef _ENDIANHELP_H
#define _ENDIANHELP_H



#include "endian_helper.h"
#include "stdint.h"

namespace endian{


    //this file contains definitions for endian functionality
    //NOTE: This only fixes endians, and is not complete for all
    //architectures. It'll work where short is 2 bytes and int is 4 bytes
    //If that's not the case, it'll probably explode.

    //a BIG thanks to the Quake 2 source code here
    extern bool BigEndianSystem;
    void InitEndian( void );		//makes use of a clever trick in Quake 2

    uint16_t ShortSwap( uint16_t s );
    uint16_t ShortNoSwap( uint16_t s );

    uint32_t LongSwap( uint32_t i );
    uint32_t LongNoSwap( uint32_t i );

    uint64_t Long64Swap( uint64_t i );
    uint64_t LongNo64Swap( uint64_t i );


    float FloatSwap( float f );
    float FloatNoSwap( float f );

    int8_t ByteSwap( int8_t c );
    int8_t ByteNoSwap( int8_t c );

    //Use these functions
    extern uint16_t (*BigShort) ( uint16_t s );
    extern uint16_t (*LittleShort) ( uint16_t s );
    extern uint32_t (*BigLong) ( uint32_t i );
    extern uint32_t (*LittleLong) ( uint32_t i );
    extern uint64_t (*BigDouble) ( uint64_t i );
    extern uint64_t (*LittleDouble) ( uint64_t i );

    extern int8_t (*LittleByte) ( int8_t c );
    extern int8_t (*BigByte) ( int8_t c );


    extern float (*BigFloat) ( float f );
    extern float (*LittleFloat) ( float f );

}

#endif
