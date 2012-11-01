#include "uParserSam.h"
#include "../utility/utility.h"
#include <sstream>
using namespace std;

namespace NGS
{

uParserSam::uParserSam(): uParserBase()
{

}

uParserSam::~uParserSam()
{

}

void uParserSam::init(const std::string& filename, bool header )
{
    uParserBase::init(filename, header);
    _parseHeader();
}

void uParserSam::init(std::iostream* stream, bool header )
{
    uParserBase::init(stream, header);
    _parseHeader();
}

uToken uParserSam::getNextEntry()
{
    try
    {
        char line[4096];
        if (m_pIostream->getline(line, 4096))
        {
            std::stringstream ss;
            ss << line;
            //String to test for int
            std::string flag;
            std::string start_pos;
            std::string MAPQual;
            std::string pNext;
            std::string Tlen;
            std::string chr;
            std::string end_pos;
            std::string score;
            std::string seq_name;
            std::string qual;
            std::string seq;
            std::string cigar;
            std::string RNext;
            std::stringstream token_infos;

            ss >> seq_name >> flag >> chr >> start_pos >> MAPQual >> cigar>>RNext>>pNext>>Tlen>>seq>>qual;
            token_infos << "CHR\t" << chr << "\n";
            token_infos << "START_POS\t" << start_pos << "\n";
            token_infos << "FLAGS\t" << flag << "\n";
            token_infos << "SEQ_NAME\t" << seq_name << "\n";
            token_infos << "MAP_SCORE\t" << MAPQual << "\n";
            token_infos << "SEQUENCE\t" << seq << "\n";
            token_infos << "CIGAR\t" << cigar << "\n";
            token_infos << "PHRED_SCORE\t" << qual << "\n";


            std::string strand="+";
            if (utility::querySamFlag(std::stoi(flag),SamQuery::SEQ_REV_STRAND))
                strand="-";

             token_infos << "STRAND\t" << strand << "\n";


           /**< Currently ignore aligment tags */
           //TODO support aligment tags

           // if (!ss.eof())
           //     throw uParser_invalid_line()<<string_error("Invalid line in Sam file, superfluous lines final token \n");
            uToken token(token_infos);
            return token;
        }
        else
        {
#ifdef DEBUG
            std::cerr << "Reached end of file." << std::endl;
#endif
            end_of_file_throw e;
            e << string_error("Reached end of file.");
            throw e;
        }
    }
    catch(invalid_uToken_throw& e)
    {
        throw e;
    }
}


void uParserSam::_processSamline(std::stringstream & curSStream)
{

}

/** \brief Parse ou same header, if any. Store in member variables
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
        if (temp.find("@HD")!=string::npos)
        {

            string data;
            string format,sortType;
            bool VN=false, SORT=false;
            while (!Infostream.eof())
            {
                Infostream >>data;
                /**< Format Version */
                if (data.find("VN:")!=string::npos)
                {
                    if (VN)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple VN tags in @HD line: \n"+lineString);
                    VN=true;
                    data.erase(0,3);
                    format=data;
                }
                /**< Sorting order */
                else if (data.find("SO:")!=string::npos)
                {
                    if (SORT)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple SO tags in @HD line: \n"+lineString);
                    //TODO validate sort
                    SORT=true;
                    data.erase(0,3);
                    sortType=data;
                }
                else
                    throw uParser_invalid_Sam_header()<<string_error("Invalid sam header line: \n"+lineString);
            }
            if (!VN)
                throw uParser_invalid_Sam_header()<<string_error("Missing VN tag in @HD header, failling: \n"+lineString);


        }
        else if (temp.find("@RG")!=string::npos)
        {

        }
        else if (temp.find("@PG")!=string::npos)
        {

        }
        /**< Reference sequence dictionnary line */
        else if (temp.find("@SQ")!=string::npos)
        {
            string data;
            string chrom;
            string REF,AssID,MD5,Species,URI;
            long long int refSeqlenght;
            bool SN=false, LN=false, AS=false, M5=false, SP=false, UR=false;
            while (!(Infostream.eof()))
            {
                Infostream >> data;
                /**< Reference sequence name */
                if (data.find("SN:")!=string::npos)
                {
                    if (SN)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple SN tag in @SQ header, failling: \n"+lineString);
                    SN=true;
                    data.erase(0,3);
                    chrom=data;
                }
                /**< Reference sequence lenght */
                else if (data.find("LN:")!=string::npos)
                {
                    if (LN)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple LN tag in @SQ header, failling: \n"+lineString);
                    LN=true;
                    data.erase(0,3);
                    refSeqlenght=std::stoll(data);
                }
                else if (data.find("AS:")!=string::npos)
                {
                    if (AS)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple AS tag in @SQ header, failling: \n"+lineString);
                    AS=true;
                    data.erase(0,3);
                    AssID=data;
                }
                else if (data.find("M5:")!=string::npos)
                {
                    if (M5)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple M5 tag in @SQ header, failling: \n"+lineString);
                    M5=true;
                    data.erase(0,3);
                    MD5=data;
                }
                else if (data.find("SP:")!=string::npos)
                {
                    if (SP)
                        throw uParser_invalid_Sam_header()<<string_error("Multiple SP tag in @SQ header, failling: \n"+lineString);
                    SP=true;
                    data.erase(0,3);
                    Species=data;
                }
                else if (data.find("UR:")!=string::npos)
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
            m_headerData._addToParam(header_param::CHR_SIZE,std::to_string(refSeqlenght));
        } /**< Invalid, fail */
        else
        {
            throw uParser_invalid_Sam_header()<<string_error("Invalid sam header line: \n"+lineString);

        }
    }

}

DerivedParserRegister<uParserSam> uParserSam::reg("SAM");

}
