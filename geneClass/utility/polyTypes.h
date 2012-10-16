#ifndef POLYTYPES_H_INCLUDED
#define POLYTYPES_H_INCLUDED

#include <iostream>
#include <map>
#include <string>
const char* EMPTY_STRING = "";

enum ValueTypes {
	TYPE_NONE,
#define X(a,b,t) t,
	X( char, asChar, TYPE_CHAR)
	X( unsigned char, asUChar, TYPE_UCHAR )
	X( short, asShort, TYPE_SHORT )
	X( unsigned short, asUShort, TYPE_USHORT )
	X( int, asInt, TYPE_INT )
	X( unsigned int, asUInt, TYPE_UINT )
	X( long, asLong, TYPE_LONG )
	X( unsigned long, asULong, TYPE_ULONG )
	X( float, asFloat, TYPE_FLOAT )
	X( double, asDouble, TYPE_DOUBLE )
	X( char*, asStr, TYPE_STR )
#undef X
};

union Values {
#define X(a,b,t) a b; \
Values(a in) { b = in; }\
operator a() { return b; }
	X( char, asChar, TYPE_CHAR)
	X( unsigned char, asUChar, TYPE_UCHAR )
	X( short, asShort, TYPE_SHORT )
	X( unsigned short, asUShort, TYPE_USHORT )
	X( int, asInt, TYPE_INT )
	X( unsigned int, asUInt, TYPE_UINT )
	X( long, asLong, TYPE_LONG )
	X( unsigned long, asULong, TYPE_ULONG )
	X( float, asFloat, TYPE_FLOAT )
	X( double, asDouble, TYPE_DOUBLE )
	X( char*, asStr, TYPE_STR )
#undef X
	Values() { asULong = 0; }
};

class PolyType {
	Values value;
	ValueTypes type;
public:
	PolyType() : value(0L), type(TYPE_NONE) { }
#define X(a,b,t) a b() { return ((type == t) ? (a)value : (a)0); } \
PolyType(a in) { value =  Values(in); type = t; } \
operator a() { return ((type == t) ? (a)value : (a)0); }
	X( char, asChar, TYPE_CHAR)
	X( unsigned char, asUChar, TYPE_UCHAR )
	X( short, asShort, TYPE_SHORT )
	X( unsigned short, asUShort, TYPE_USHORT )
	X( int, asInt, TYPE_INT )
	X( unsigned int, asUInt, TYPE_UINT )
	X( long, asLong, TYPE_LONG )
	X( unsigned long, asULong, TYPE_ULONG )
	X( float, asFloat, TYPE_FLOAT )
	X( double, asDouble, TYPE_DOUBLE )
//	X( char*, asStr, TYPE_STR ) //special case required
#undef X
	//This special case in line with the intended use of this class...
	//  if the type is not a string then it returns a pointer to an empty string rather than
	//  a cast to zero.
	char* asStr() { return ((type == TYPE_STR) ? (char*)value : (char*)EMPTY_STRING); }
	PolyType(char* in) { value =  Values(in); type = TYPE_STR; }
	operator char*() { return ((type == TYPE_STR) ? (char*)value : (char*)EMPTY_STRING); }

	friend std::ostream& operator<<(std::ostream& out, PolyType poly);
};

std::ostream& operator<<(std::ostream& out, PolyType poly) {
	switch(poly.type) {
	#define X(a,b,t) case t: { out << poly.value.b ; break; }
		X( char, asChar, TYPE_CHAR)
		X( unsigned char, asUChar, TYPE_UCHAR )
		X( short, asShort, TYPE_SHORT )
		X( unsigned short, asUShort, TYPE_USHORT )
		X( int, asInt, TYPE_INT )
		X( unsigned int, asUInt, TYPE_UINT )
		X( long, asLong, TYPE_LONG )
		X( unsigned long, asULong, TYPE_ULONG )
		X( float, asFloat, TYPE_FLOAT )
		X( double, asDouble, TYPE_DOUBLE )
		X( char*, asStr, TYPE_STR )
	#undef X
	}
	return out;
}

#endif // POLYTYPES_H_INCLUDED
