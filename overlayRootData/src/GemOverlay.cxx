//
// The GemOverlay class contains GEM information for Overlays
//

#include "overlayRootData/GemOverlay.h"
#include "Riostream.h"

ClassImp(GemOverlayTileList)
ClassImp(GemOverlayOnePpsTime)
ClassImp(GemOverlayCondArrivalTime)
ClassImp(GemOverlay)

void GemOverlayTileList::Print(Option_t *option) const
{
    using namespace std;
    TObject::Print(option);
}

void GemOverlayOnePpsTime::Print(Option_t *option) const 
{
    using namespace std;
    TObject::Print(option);
}

GemOverlay::GemOverlay() 
{
    Clear();
}

GemOverlay::GemOverlay(const GemOverlay& GemOverlay) :TObject(GemOverlay) 
{
    m_tkrVector         = GemOverlay.m_tkrVector;
    m_roiVector         = GemOverlay.m_roiVector;
    m_cal_Le_Vector     = GemOverlay.m_cal_Le_Vector;
    m_cal_He_Vector     = GemOverlay.m_cal_He_Vector;
    m_cno_Vector        = GemOverlay.m_cno_Vector;
    m_conditionSummary  = GemOverlay.m_conditionSummary;
    m_tileList          = GemOverlay.m_tileList;
    m_liveTime          = GemOverlay.m_liveTime;
    m_prescaled         = GemOverlay.m_prescaled;
    m_sent              = GemOverlay.m_sent;
    m_discarded         = GemOverlay.m_discarded;
    m_onePpsTime        = GemOverlay.m_onePpsTime;
    m_triggerTime       = GemOverlay.m_triggerTime;
    m_deltaEventTime    = GemOverlay.m_deltaEventTime;
    m_deltaWindOpenTime = GemOverlay.m_deltaWindOpenTime;

    m_condArr.init(GemOverlay.m_condArr.condArr());
    m_missed = GemOverlay.m_missed;
}

void GemOverlay::initTrigger(UShort_t tkr, UShort_t roi,
                     UShort_t calLE,
                     UShort_t calHE, UShort_t cno,
                     UShort_t conditionSummary,
                     UShort_t missed,
                     const GemOverlayTileList &tileList) 
{
    m_tkrVector        = tkr;
    m_roiVector        = roi;
    m_cal_Le_Vector    = calLE;
    m_cal_He_Vector    = calHE;
    m_cno_Vector       = cno;
    m_conditionSummary = conditionSummary;
    m_missed           = missed;
    m_tileList         = tileList;
}

void GemOverlay::initSummary(UInt_t liveTime, UInt_t prescaled,
                     UInt_t discarded, UInt_t condArrTime,
                     UInt_t triggerTime, const GemOverlayOnePpsTime &time,
                     UShort_t deltaEventTime, UShort_t deltaWindOpenTime)
{
    m_liveTime          = liveTime;
    m_prescaled         = prescaled;
    m_discarded         = discarded;
    m_sent              = 0;
    m_triggerTime       = triggerTime;
    m_onePpsTime        = time;
    m_deltaEventTime    = deltaEventTime;
    m_deltaWindOpenTime = deltaWindOpenTime;

    m_condArr.init(condArrTime);
}

void GemOverlay::Clear(Option_t *option) 
{
    m_tkrVector         = 0;   
    m_roiVector         = 0; 
    m_cal_Le_Vector     = 0;
    m_cal_He_Vector     = 0;
    m_cno_Vector        = 0;
    m_conditionSummary  = 0;
    m_liveTime          = 0;   
    m_prescaled         = 0;
    m_discarded         = 0;
    m_sent              = 0;
    m_triggerTime       = 0;
    m_deltaEventTime    = 0;
    m_deltaWindOpenTime = 0;
    m_missed            = 0;

    m_tileList.Clear(option);
    m_condArr.init(0);
    m_onePpsTime.Clear(option);
}

void GemOverlay::Print(Option_t *option) const 
{
    using namespace std;
    TObject::Print(option);
    //cout.precision(2);
}
