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
/**< Modified by Alexei Nordell Markovits, August
Uses fixed bit sizes to faciliate reading from fixed size binary file.
Added double overload.
Added char Overload
 */

#include "endian_helper.h"
#include <iostream>
namespace endian{

typedef unsigned char byte;

bool BigEndianSystem;

uint16_t(*BigShort) ( uint16_t s );
uint16_t(*LittleShort) ( uint16_t s );
uint32_t(*BigLong) ( uint32_t i );
uint32_t(*LittleLong) ( uint32_t i );
uint64_t(*BigDouble) ( uint64_t i );
uint64_t(*LittleDouble) ( uint64_t i );


float (*BigFloat) ( float f );
float (*LittleFloat) ( float f );


int8_t (*BigByte) ( int8_t c );
int8_t (*LittleByte) ( int8_t c );


int8_t ByteSwap( int8_t s )
{
//	byte b1;

	//b1 = s & 255;
	//b2 = (s >> 8) & 255;

	//return (b1 << 8) + b2;
	return s;
}

int8_t ByteNoSwap( int8_t c )
{
	return c;
}

//adapted from Quake 2 source

uint16_t ShortSwap( uint16_t s )
{
	byte b1, b2;

	b1 = s & 255;
	b2 = (s >> 8) & 255;

	return ((uint16_t)b1 << 8) + b2;
}

uint16_t ShortNoSwap( uint16_t s )
{
	return s;
}



uint32_t LongSwap (uint32_t i)
{
	byte b1, b2, b3, b4;
	b1 = i & 255;
	b2 = ( i >> 8 ) & 255;
	b3 = ( i>>16 ) & 255;
	b4 = ( i>>24 ) & 255;

	return ((uint32_t)b1 << 24) + ((uint32_t)b2 << 16) + ((uint32_t)b3 << 8) + b4;
}

uint32_t LongNoSwap( uint32_t i )
{

	return i;
}


uint64_t Long64Swap (uint64_t i)
{
	byte b1, b2, b3, b4,b5,b6,b7,b8;
	b1 = i & 255;
	b2 = ( i >> 8 ) & 255;
	b3 = ( i>>16 ) & 255;
	b4 = ( i>>24 ) & 255;
	b5 = ( i>>32 ) & 255;
	b6 = ( i>>40 ) & 255;
	b7 = ( i>>48 ) & 255;
	b8 = ( i>>56 ) & 255;


	return ((uint64_t)b1 << 56) + ((uint64_t)b2 << 48) + ((uint64_t)b3 << 40) +((uint64_t)b4 << 32) +((uint64_t)b5 << 24) +((uint64_t)b6 << 16) + ((uint64_t)b7 << 8)+ b8;
}

uint64_t LongNo64Swap( uint64_t i )
{
	return i;
}




float FloatSwap( float f )
{
	union
	{
		float f;
		byte b[4];
	} dat1, dat2;

	dat1.f = f;
	dat2.b[0] = dat1.b[3];
	dat2.b[1] = dat1.b[2];
	dat2.b[2] = dat1.b[1];
	dat2.b[3] = dat1.b[0];
	return dat2.f;
}

float FloatNoSwap( float f )
{

	return f;
}


void InitEndian( void )
{
	//clever little trick from Quake 2 to determine the endian
	//of the current system without depending on a preprocessor define

	byte SwapTest[2] = { 1, 0 };

	if( *(short *) SwapTest == 1 )
	{
		//little endian
		BigEndianSystem = false;

		//set func pointers to correct funcs
		BigShort = ShortSwap;
		LittleShort = ShortNoSwap;
		BigLong = LongSwap;
		LittleLong = LongNoSwap;
		BigFloat = FloatSwap;
		LittleFloat = FloatNoSwap;
        BigDouble = Long64Swap;
		LittleDouble = LongNo64Swap;

        BigFloat = FloatSwap;
		LittleFloat = FloatNoSwap;

        BigByte = ByteSwap;
        LittleByte = ByteNoSwap;

	}
	else
	{
		//big endian
		BigEndianSystem = true;

		BigShort = ShortNoSwap;
		LittleShort = ShortSwap;
		BigLong = LongNoSwap;
		LittleLong = LongSwap;
		BigDouble = LongNo64Swap;
		LittleDouble = Long64Swap;

		BigFloat = FloatNoSwap;
		LittleFloat = FloatSwap;

        BigByte = ByteNoSwap;
        LittleByte = ByteSwap;

	}
}

}

