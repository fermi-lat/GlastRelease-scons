// $Header$
#include "IRFConverter.h"


#include "instrument/Scintillator.h"
#include "instrument/SiDetector.h"
#include "instrument/CsIDetector.h"
#include "instrument/MCTruth.h"
#include "instrument/Tower.h"

// note initialization of number of layers
IRFConverter::IRFConverter():m_cal(8),m_tkr(18){}

void IRFConverter::forward (const Tower& t) {
    m_id = t.getId();
}
void IRFConverter::forward ( const CsIDetector& csi_det) {
    if ( !(csi_det.empty()) ) {
        m_cal.load(csi_det, m_id);
    }
}
void IRFConverter::forward ( const SiDetector& si_det) {
    if ( !(si_det.empty()) ) {
        m_tkr.load(si_det, m_id);
    }
}
void IRFConverter::forward ( const Scintillator& scint) {
    if ( !(scint.empty()) ) {
        m_acd.load(scint);
    }
}
void IRFConverter::forward ( const MCTruth& mc){
}

void IRFConverter::printOn(std::ostream& out) const
{
   m_acd.printOn(out);
   m_tkr.printOn(out);
   m_cal.printOn(out);
}

void IRFConverter::clear()
{
    m_acd.clear();
    m_tkr.clear();
    m_cal.clear();
}
