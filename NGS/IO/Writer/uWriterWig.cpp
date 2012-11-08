#include "uWriterWig.h"

namespace NGS {
using namespace std;
/** \brief Print the values of a token in SAM format in current file
  * \param const uToken& token: the token to print.
   */
  void uWriterWig::writeToken(const uToken& token) {
 	try {
	//	*m_pOstream <<seq_name <<TAB << flag <<TAB << chr <<TAB <<start_pos<<TAB <<MAPQual<<TAB <<cigar<<TAB << RNext<<TAB <<pNext<<TAB <<TLEN<<TAB<<seq<<TAB<<PHREDQual<<endl;
	}
	catch(param_not_found& e) {
		throw e;
	}
	*m_pOstream << std::endl;
}

DerivedRegister<uWriterWig> uWriterWig::reg("WIG");
} // End of namespace NGS
