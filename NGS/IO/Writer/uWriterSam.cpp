
#include "uWriterSam.h"

namespace NGS {
using namespace std;

/** \brief Queries our header object and if the appropriate information is there, write the reference sizes
 *  Also output a mostly meaningless @HQ h
 * \return void
 *
 */
void uWriterSam::writeHeader(){
    const string TAB="\t";
    const string SN="SN:";
    const string SQ="LN:";
    /**< Write Version number and sort guarantee */
    *m_pOstream<<"@HD    VN:1.0  SO:unknown"<<std::endl;
  if (m_headerData.isParamSet(header_param::CHR))
    {
        auto chrVector=m_headerData.getParamVector(header_param::CHR);
        auto chrLenghts= m_headerData.getParamVector(header_param::CHR_SIZE);

        /**< Write every header */
        for (unsigned int i=0; i<chrVector.size();i++)
            {
                *m_pOstream<<"@SQ"<<TAB<<SN<<chrVector.at(i)<<TAB<<SQ<<chrLenghts.at(i)<<"\n";
            }
    }
}

/** \brief Print the values of a token in SAM format in current file
  * \param const uToken& token: the token to print.
  * \brief Note, we allow an element to map outside of the reference, has there is some contradiction in how tools manage this.
  */
void uWriterSam::writeToken(const uToken& token) {
	try {
        const string TAB="\t";
		/**< If score or strand are not set, we replace them with the value "." */
        std::string flag="0";
        std::string start_pos="0";
        std::string MAPQual="0";
        std::string pNext="*";
        std::string Tlen="0";
        std::string chr="*";
        std::string end_pos;
        std::string PHREDQual="*";
        std::string seq_name="*";
        std::string qual="0";
        std::string seq="*";
        std::string cigar="*";
        std::string RNext="*";
        std::string TLEN="0";
		if (token.isParamSet(token_param::CHR)) {
			chr = token.getParam(token_param::CHR);
		}
		if(chrValidate(chr)==false)
            throw uWriter_missing_mandatory_header()<<string_error("Sam entry has no-corresponding RNAME");

        if (token.isParamSet(token_param::START_POS)) {
			start_pos = token.getParam(token_param::START_POS);
		}
        if (token.isParamSet(token_param::END_POS)) {
			end_pos = token.getParam(token_param::END_POS);
		}
        if (token.isParamSet(token_param::FLAGS)) {
			flag = token.getParam(token_param::FLAGS);
		}
        if (token.isParamSet(token_param::SEQUENCE)) {
			seq = token.getParam(token_param::SEQUENCE);
		}
        if (token.isParamSet(token_param::CIGAR)) {
			cigar = token.getParam(token_param::CIGAR);
		}
        if (token.isParamSet(token_param::PHRED_SCORE)) {
			PHREDQual = token.getParam(token_param::PHRED_SCORE);
		}
        if (token.isParamSet(token_param::SEQ_NAME)) {
			seq_name = token.getParam(token_param::SEQ_NAME);
		}
        if (token.isParamSet(token_param::MAP_SCORE)) {
			MAPQual = token.getParam(token_param::MAP_SCORE);
		}
        if (token.isParamSet(token_param::FLAGS)) {
			flag = token.getParam(token_param::FLAGS);
		}
         /**< End pos but no way to get that size, infer CIGAR */
        if ((cigar=="*")&&(seq=="*")&&(end_pos!="0"))
        {
            cigar=utility::to_string ( (utility::stoi(end_pos)-utility::stoi(start_pos) )+1 )+"M";
        }
		*m_pOstream <<seq_name <<TAB << flag <<TAB << chr <<TAB <<start_pos<<TAB <<MAPQual<<TAB <<cigar<<TAB << RNext<<TAB <<pNext<<TAB <<TLEN<<TAB<<seq<<TAB<<PHREDQual<<"\n";
	}
	catch(param_not_found& e) {
		throw e;
	}
    catch(...) {
		throw ;
	}
}
 /** \brief True if found in header list or no header list
  *
  * \param chr const std::string& chromose to check for
  * \return bool Valid or not
  *
  */
 bool uWriterSam::chrValidate(const std::string& chr){
    if (m_headerData.isParamSet(header_param::CHR))
    {
        auto chrVector=m_headerData.getParamVector(header_param::CHR);
        if (std::find(chrVector.begin(), chrVector.end(), chr)!=chrVector.end() )
            return true;
        return false;
    }
return true;
}
//DerivedRegister<uWriterSam> uWriterSam::reg("SAM");
} // End of namespace NGS
