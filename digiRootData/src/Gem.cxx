//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//

#include "digiRootData/Gem.h"
#include <iostream>

ClassImp(GemTileList)
ClassImp(GemOnePpsTime)
ClassImp(Gem)

void GemTileList::Print(Option_t *option) const{
    using namespace std;
    TObject::Print(option);

}

void GemOnePpsTime::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);

}

Gem::Gem() {
    Clear();
}

Gem::Gem(const Gem& gem) :TObject(gem) {

    m_tkrVector = gem.m_tkrVector;
    m_roiVector = gem.m_roiVector;
    m_cal_Le_Vector = gem.m_cal_Le_Vector;
    m_cal_He_Vector = gem.m_cal_He_Vector;
    m_cno_Vector = gem.m_cno_Vector;
    m_conditionSummary = gem.m_conditionSummary;
    m_tileList = gem.m_tileList;
    m_liveTime = gem.m_liveTime;
    m_prescaled = gem.m_prescaled;
    m_sent = gem.m_sent;
    m_discarded = gem.m_discarded;
    m_onePpsTime = gem.m_onePpsTime;
    m_triggerTime = gem.m_triggerTime;
    m_deltaEventTime = gem.m_deltaEventTime;
}

void Gem::initTrigger(UShort_t tkr, UShort_t roi,
                     UShort_t calLE,
                     UShort_t calHE, UShort_t cno,
                     UShort_t conditionSummary,
                     const GemTileList &tileList) 
{
    m_tkrVector = tkr;
    m_roiVector = roi;
    m_cal_Le_Vector = calLE;
    m_cal_He_Vector = calHE;
    m_cno_Vector = cno;
    m_conditionSummary = conditionSummary;
    m_tileList = tileList;
}

void Gem::initSummary(UInt_t liveTime, UInt_t prescaled,
                     UInt_t discarded, UInt_t sent,
                     UInt_t triggerTime, const GemOnePpsTime &time,
                     UInt_t deltaEventTime)
{
    m_liveTime = liveTime;
    m_prescaled = prescaled;
    m_discarded = discarded;
    m_sent = sent;
    m_triggerTime = triggerTime;
    m_onePpsTime = time;
    m_deltaEventTime = deltaEventTime;

}

void Gem::Clear(Option_t *option) {

    m_tkrVector = 0;   
    m_roiVector = 0; 
    m_cal_Le_Vector = 0;
    m_cal_He_Vector = 0;
    m_cno_Vector = 0;
    m_conditionSummary = 0;
    m_tileList.Clear(option);
    m_liveTime = 0;   
    m_prescaled = 0;
    m_discarded = 0;
    m_sent = 0;
    m_triggerTime = 0;
    m_onePpsTime.Clear(option);
    m_deltaEventTime = 0;
}

void Gem::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    //cout.precision(2);
}
