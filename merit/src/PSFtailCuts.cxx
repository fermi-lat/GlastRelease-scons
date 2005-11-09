// $Header$
// PSFtailCuts.cxx: implementation of the PSFtailCuts class.
//
//////////////////////////////////////////////////////////////////////

#include "PSFtailCuts.h"
#include "Cut.h"
#include <sstream>
#include <cmath>
//=============================================================================
class AbsValueCut : public Analyze {
    friend class PSFtailCuts;
    AbsValueCut(const Tuple& t, const std::string& name, double value)
        : Analyze(t,name)
        , m_value(value)
    {
        std::stringstream label;
        label << "abs(" << name << ")<" << value;
        set_name(label.str());
    }

    virtual bool apply() {
        return fabs(item()) < m_value;
    }
    double m_value;
};


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PSFtailCuts::PSFtailCuts(const Tuple& t): AnalysisList("CT PSF tail cuts")
{
    try {
    push_back( new Cut(t, "CTgoodCal>0.25") );
    }catch(...){}

    push_back( new Cut(t, "Tkr1Zdir<-0.25") );
}
