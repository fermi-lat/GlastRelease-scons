#include "digiRootData/LiveTime.h"

ClassImp(LiveTime)

LiveTime::LiveTime() {
    // Default constructor.
    // Assign default values to members
    Clean();
}

LiveTime::~LiveTime() {
  // Destructor    
}

void LiveTime::Clean(Option_t *option) {
    m_l1tCount = 0;  
    m_ppsCount = 0;  

    m_ppsTime = 0;
    m_ppsTimeOverFlow = 0;
    
    m_liveTime = 0;
    m_liveTimeOverFlow = 0;

    m_runTime = 0;
    m_runTimeOverFlow = 0;

    m_runCount = 0;
}


