#include "uWriterGFF.h"
namespace NGS
{
/** \brief Print the values of a token in GFF format in current file
  * \param const uToken& token: the token to print.
  */


/* Note, there are ambiguities in the GFF format. The spec at http://www.sanger.ac.uk/resources/software/gff/spec.html
describes the first column as

"The name of the sequence. Having an explicit sequence name allows a feature file to be prepared for a data set of multiple sequences.
Normally the seqname will be the identifier of the sequence in an accompanying fasta format file. An alternative is that <seqname>
is the identifier for a sequence in a public database, such as an EMBL/Genbank/DDBJ accession number. Which is the case, and which
file or database to use, should be explained in accompanying information. "


However, practically, use tend to set column one as scaffp;d ID ( ex: chromosome ). This same langage is used for GFF3

As such, we write SEQ_NAME in the first column, if available. HOwever, if unavailable, we write CHR in first column.

Note that the GFF3 parsers set CHR as first name always.
 */
void uWriterGFF::writeToken(const uToken& token)
{
    try
    {
        /**< If score or strand are not set, we replace them with the value "." */
        const std::string TAB="\t";
        std::string seqname=".";
        std::string feature=".";
        std::string score=".";
        std::string strand=".";
        std::string start;
		std::string end;
		std::string source=".";
		std::string phase=".";
		std::string extra;

		/**< No default values */
		start=token.getParam(token_param::START_POS);
		end=token.getParam(token_param::END_POS);

        if (token.isParamSet(token_param::SCORE))
       	 { score = token.getParam(token_param::SCORE); }
       	 /**< Practically speaking, most people seem to use the chr as the seqname of a GFF. So if one is not evailable, we use  chr */
		if (token.isParamSet(token_param::SEQ_NAME))
       	 { seqname = token.getParam(token_param::SEQ_NAME); }
		else
		if(token.isParamSet(token_param::CHR))
       		 { seqname = token.getParam(token_param::CHR); }

		if(token.isParamSet(token_param::FEATURE_TYPE))
			{ feature = token.getParam(token_param::FEATURE_TYPE); }
        /**< If score is not set, but strand is, we replace score by "." */
        if (token.isParamSet(token_param::STRAND))
        {  strand = token.getParam(token_param::STRAND);   }

		if (token.isParamSet(token_param::SOURCE))
        {  source = token.getParam(token_param::SOURCE);   }

		if(token.isParamSet(token_param::PHASE))
			{ phase = token.getParam(token_param::PHASE); }

		if(token.isParamSet(token_param::EXTRA))
       		 { extra = token.getParam(token_param::EXTRA); }

        *m_pOstream << seqname<<TAB<< source<<TAB << feature<<TAB << start<<TAB<<end <<TAB << score<<TAB<<strand <<TAB<< phase<<std::endl;
    }
    catch(param_not_found& e)
    {
        throw e;
    }
}
//DerivedRegister<uWriterGFF> uWriterGFF::reg("GFF");
} // End of namespace NGS
