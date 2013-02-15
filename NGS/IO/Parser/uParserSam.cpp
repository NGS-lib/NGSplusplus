#include "uParserSam.h"
#include "../utility/utility.h"
#include <sstream>
//#include "uParserFactory.h"
namespace NGS
{

using namespace boost::xpressive;
/** \brief Default empty constructor
 */
uParserSam::uParserSam()
:uParserBase()
,s10(10)
,s11(11)
,s12(12)
{

}
/** \brief Default empty destructor
 */
uParserSam::~uParserSam()
{
}

/** \brief Called when created, loads stream and parses Sam file header, loading mandatory info
 *
 * \param filename const std::string& Path to file
 * \param header bool leave at false, ignored if set
 * \return void
 *
 */
void uParserSam::init(const std::string& filename, bool header )
{
    uParserBase::init(filename, header);
    _parseHeader();


}

/** \brief Called when created, reads from stream and parses Sam file header, loading mandatory info
 *
 * \param stream std::istream*
 * \param header bool leave at false, ignored if set
 * \return void
 *
 */
void uParserSam::init(std::istream* stream, bool header )
{
    uParserBase::init(stream, header);
    _parseHeader();

}

/**< Old parsing code, preserved for now */

///** \brief Parses next line in Sam file
// *
// * \return uToken Token created from Sam
// *
// */
//uToken uParserSam::getNextEntryWithRegex()
//{
//    try
//    {
//        char line[4096];
//        if (m_pIostream->getline(line, 4096))
//        {
//            std::stringstream ss;
//            ss << line;
//            //String to test for int
//            std::string flag;
//            std::string start_pos;
//            std::string MAPQual;
//            std::string pNext;
//            std::string Tlen;
//            std::string chr;
//            std::string end_pos;
//            std::string score;
//            std::string seq_name;
//            std::string qual;
//            std::string seq;
//            std::string cigar;
//            std::string RNext;
//            std::stringstream token_infos;
//            //TODO : Regex!
//            ss >> seq_name >> flag >> chr >> start_pos >> MAPQual >> cigar>>RNext>>pNext>>Tlen>>seq>>qual;
//            token_infos << "SEQ_NAME\t" << seq_name << "\n";
//            token_infos << "FLAGS\t" << flag << "\n";
//            token_infos << "CHR\t" << chr << "\n";
//            token_infos << "START_POS\t" << start_pos << "\n";
//            token_infos << "MAP_SCORE\t" << MAPQual << "\n";
//            token_infos << "SEQUENCE\t" << seq << "\n";
//            token_infos << "CIGAR\t" << cigar << "\n";
//            token_infos << "PHRED_SCORE\t" << qual << "\n";
//
//
//            std::string strand="+";
//            if (utility::querySamFlag(utility::stoi(flag),SamQuery::SEQ_REV_STRAND))
//                strand="-";
//
//            token_infos << "STRAND\t" << strand << "\n";
//
//
//            /**< Currently ignore aligment tags */
//            //TODO support aligment tags
//
//            // if (!ss.eof())
//            //     throw uParser_invalid_line()<<string_error("Invalid line in Sam file, superfluous lines final token \n");
//
//            return uToken(token_infos);
//        }
//        else
//        {
//#ifdef DEBUG
//            std::cerr << "Reached end of file." << std::endl;
//#endif
//            end_of_file_throw e;
//            e << string_error("Reached end of file.");
//            throw e;
//        }
//    }
//    catch(invalid_uToken_throw& e)
//    {
//        throw e;
//    }
//}


/** \brief Parses next line in Sam file with xpressive regex
 *
 * \return uToken Token created from Sam
 *
 */
uToken uParserSam::getNextEntry()
{
    try
    {
        std::string strLine;
        std::getline(*m_pIostream, strLine);
         m_rawString=strLine;
		//std::stringstream token_infos;
		/**< For readibility sake, macro or split this up. */
		uToken ourToken;
        if( regex_match( strLine, what, staticSam ) )
        {
			ourToken._setParamNoValidate(token_param::SEQ_NAME, what[1]);
			ourToken._setParamNoValidate(token_param::FLAGS, what[2]);
			ourToken._setParamNoValidate(token_param::CHR, what[3]);
			ourToken._setParamNoValidate(token_param::START_POS, what[4]);
			ourToken._setParamNoValidate(token_param::MAP_SCORE, what[5]);
			ourToken._setParamNoValidate(token_param::CIGAR, what[6]);
			ourToken._setParamNoValidate(token_param::TEMPLATE_LENGHT, what[9]);
			ourToken._setParamNoValidate(token_param::SEQUENCE, what[10]);
			ourToken._setParamNoValidate(token_param::PHRED_SCORE, what[11]);
       //     token_infos << "SEQ_NAME\t" << what[1] << "\n";
       //     token_infos << "FLAGS\t" <<  what[2] << "\n";
       //     token_infos << "CHR\t" << what[3] << "\n";
		//	token_infos << "START_POS\t" << what[4] << "\n";
		//	token_infos << "MAP_SCORE\t" << what[5] << "\n";
		//	token_infos << "CIGAR\t" << what[6] << "\n";
		//	token_infos << "TEMPLATE_LENGHT\t" << what[9] << "\n";
			/**< Skip RNEXT and PNEXT and Template Lenght */
      //      token_infos << "SEQUENCE\t" << what[10] << "\n";
      //    	token_infos << "PHRED_SCORE\t" << what[11] << "\n";
        }
        else
        {
            throw uParser_invalid_Sam_line()<<string_error("SAM line, failling validation. Line is:\n"+strLine);
        }
        std::string strand="+";
        if (utility::SAM::querySamFlag(utility::stoi(what[2]),SamQuery::SEQ_REV_STRAND))
            strand="-";
		ourToken._setParamNoValidate(token_param::STRAND, strand);
		//token_infos << "STRAND\t" << strand << "\n";
		return ourToken;
		//return uToken(token_infos,false,false);
        }
catch(invalid_uToken_throw& e)
{
    throw e;
}
}




/** \brief Parse SAM header, if any. Store in member variables
 *
 * \return void
 *
 */
void uParserSam::_parseHeader()
{
    std::string lineString;
    /**< If header, otherwise skip */
    while (m_pIostream->peek()=='@')
    {
        getline(*m_pIostream, lineString);

        /**< Parse the line */
        std::stringstream Infostream;
        std::string temp,chrom, size;
        Infostream.str(lineString);
        /**< get Token */
        Infostream >>temp;

        /**< Header line */
        if (temp.find("@HD")!=std::string::npos)
        {

            std::string data;
            std::string format,sortType;
            bool VN=false, SORT=false;
            while (!Infostream.eof())
            {
                Infostream >>data;
                /**< Format Version */
                if (data.find("VN:")!=std::string::npos)
                {
                    if (VN)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple VN tags in @HD line: \n"+lineString);
                    VN=true;
                    data.erase(0,3);
                    format=data;
                }
                /**< Sorting order */
                else if (data.find("SO:")!=std::string::npos)
                {
                    if (SORT)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple SO tags in @SO line: \n"+lineString);
                    SORT=true;
                    data.erase(0,3);

                    if ((data!="unsorted")&&(data!="queryname")&&(data!="unknown")&&(data!="coordinate"))
                         throw uParser_invalid_Sam_header()<<string_error("Invalid sorting value in @SO line: \n"+lineString);
                    sortType=data;
                }
                else
                    throw uParser_invalid_Sam_header()<<string_error("Invalid sam header line: \n"+lineString);
            }
            if (!VN)
                throw uParser_invalid_Sam_header()<<string_error("Missing VN tag in @HD header, failling: \n"+lineString);


        }
        else if (temp.find("@RG")!=std::string::npos)
        {

        }
        else if (temp.find("@PG")!=std::string::npos)
        {

        }
        /**< Reference sequence dictionnary line */
        else if (temp.find("@SQ")!=std::string::npos)
        {
            std::string data;
            std::string chrom;
            std::string REF,AssID,MD5,Species,URI;
            long long int refSeqlenght;
            bool SN=false, LN=false, AS=false, M5=false, SP=false, UR=false;
            while (!(Infostream.eof()))
            {
                Infostream >> data;
                /**< Reference sequence name */
                if (data.find("SN:")!=std::string::npos)
                {
                    if (SN)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple SN tag in @SQ header, failling: \n"+lineString);
                    SN=true;
                    data.erase(0,3);
                    chrom=data;
                }
                /**< Reference sequence lenght */
                else if (data.find("LN:")!=std::string::npos)
                {
                    if (LN)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple LN tag in @SQ header, failling: \n"+lineString);
                    LN=true;
                    data.erase(0,3);
                    refSeqlenght=utility::stoll(data);
                }
                else if (data.find("AS:")!=std::string::npos)
                {
                    if (AS)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple AS tag in @SQ header, failling: \n"+lineString);
                    AS=true;
                    data.erase(0,3);
                    AssID=data;
                }
                else if (data.find("M5:")!=std::string::npos)
                {
                    if (M5)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple M5 tag in @SQ header, failling: \n"+lineString);
                    M5=true;
                    data.erase(0,3);
                    MD5=data;
                }
                else if (data.find("SP:")!=std::string::npos)
                {
                    if (SP)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple SP tag in @SQ header, failling: \n"+lineString);
                    SP=true;
                    data.erase(0,3);
                    Species=data;
                }
                else if (data.find("UR:")!=std::string::npos)
                {
                    if (UR)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple UR tag in @SQ header, failling: \n"+lineString);
                    UR=true;
                    data.erase(0,3);
                    Species=URI;
                }
            }
            /**< Set the data */
            if((!SN)||(!LN))
                throw uParser_invalid_Sam_header()<<string_error("Missing SN or LN tag in @SQ header, failling: \n"+lineString);
            /**< Load our data */
            m_headerData._addToParam(header_param::CHR,chrom);
            m_headerData._addToParam(header_param::CHR_SIZE,utility::to_string(refSeqlenght));
        } /**< Invalid, fail */
        else
        {
            throw uParser_invalid_Sam_header()<<string_error("Invalid sam header line: \n"+lineString);

        }
    }

}
//DerivedParserRegister<uParserSam> uParserSam::reg("SAM");
}
