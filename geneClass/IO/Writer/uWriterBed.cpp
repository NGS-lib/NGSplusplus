#include "uWriterBed.h"

namespace NGS {

void uWriterBed::init(const std::string& filename) {
}

void uWriterBed::init(std::ofstream* os) {
}

void uWriterBed::writeToken(const uToken& token) {
}

DerivedRegister<uWriterBed> uWriterBed::reg("BED");
} // End of namespace NGS
