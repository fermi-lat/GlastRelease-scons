// LdCsIData.cpp: implementation of the LdCsIData class.
//
//////////////////////////////////////////////////////////////////////

#include "LdCsIData.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

LdCsIData::LdCsIData(int numLay):TdCsIData(numLay){
    ;
};

/*! Takes in a TdCsIData object makes an xtal object and pushes in back onto
    the  calorList of xtal objects. Currently called from IRFConverter
*/
void LdCsIData::load (const CsIDetector& xtal, idents::ModuleId tower)
{
    
    // The moduleId was taken out and replaced with a zero for now    
    if (xtal.hit()) {


        calorList[xtal.layer()]->push_back(Xtal(xtal.xtalPos(), xtal.energy(),tower,
            xtal.xtalId(), xtal.Lresp(), xtal.Rresp(),
				xtal.getDiodesEnergy()));

        // put the "raw" information into the basic list
        idents::XtalId id(tower, xtal.layer(), xtal.xtalId());
        int	left =  static_cast<int>(xtal.Lresp()*1e6),
            right = static_cast<int>(xtal.Rresp()*1e6);
            m_xtals[id]=std::pair<int,int>(left,right);
    }


}



