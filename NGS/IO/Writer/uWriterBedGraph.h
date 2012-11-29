#ifndef UWRITERBEDGRAPH_H_INCLUDED
#define UWRITERBEDGRAPH_H_INCLUDED

#include <iostream>
#include <fstream>
#include "uWriterBase.h"
namespace NGS {

class uWriterBedGraph : public uWriterBase {
public:
	/** \brief Empty constructor (call object with init through factory instead)
	  */
	uWriterBedGraph() {}
	void writeToken(const uToken& token);
	void writeHeader();
private:
	static DerivedRegister<uWriterBedGraph> reg;
    const std::string m_bedGraphHeader="track type=bedGraph";
}; // End of class uWriterBed4

} // End of namespace NGS
#endif // uWriterBedGraph_H_INCLUDED

