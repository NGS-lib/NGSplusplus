#include "uBasicNGS.h"
#include "uTags.h"
#include "uRegion.h"
namespace NGS {

    uBasicNGS::uBasicNGS(uTags p_tags):uGenericNGS(p_tags.getChr(),p_tags.getStart(), p_tags.getEnd(), p_tags.getStrand()){
    }

    uBasicNGS::uBasicNGS(uRegion p_region):uGenericNGS(p_region.getChr(),p_region.getStart(), p_region.getEnd(), p_region.getStrand()){
    }

    /** \brief Copy constructor
     *
     * \param copyCop : The object to instaciate  from.
     */
    /** \brief Copy constructor
     *
     * \param copyCop : The object to instaciate  from.
     */
    uBasicNGSChrom::uBasicNGSChrom(const uGenericNGSChrom<uBasicNGSChrom,uBasicNGS> & copyCop)
    {
        VecSites=copyCop.returnVecData();
        chr= copyCop.getChr();
        m_isSorted=copyCop.getSortedStatus();
        sortGetStart=copyCop.getStartFunct();
        sortGetEnd=copyCop.getEndFunct();
        m_comptFunc=copyCop.getCompFunct();
        chromSize=copyCop.getChromSize();
    }

    uBasicNGSChrom::uBasicNGSChrom(const uBasicNGSChrom& initFrom){

        VecSites=initFrom.returnVecData();
        chr= initFrom.getChr();
        m_isSorted=initFrom.m_isSorted;
        sortGetStart=initFrom.sortGetStart;
        sortGetEnd=initFrom.sortGetEnd;
        m_comptFunc=initFrom.m_comptFunc;
        chromSize=initFrom.chromSize;

    }
    uBasicNGSChrom& uBasicNGSChrom::operator=(const uBasicNGSChrom& copFrom)
    {

        VecSites=copFrom.returnVecData();
        chr= copFrom.getChr();
        m_isSorted=copFrom.m_isSorted;
        sortGetStart=copFrom.sortGetStart;
        sortGetEnd=copFrom.sortGetEnd;
        m_comptFunc=copFrom.m_comptFunc;
        chromSize=copFrom.chromSize;

        return *this;
    }


}
