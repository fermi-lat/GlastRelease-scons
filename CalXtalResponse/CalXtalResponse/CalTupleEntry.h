// LOCAL INCLUDES

// GLAST INCLUDES

// EXTLIB INCLUDES
#include "Rtypes.h"

// STD INCLUDES
#include <cstring>

/// Single entry in CalTuple
struct CalTupleEntry {
  void Clear() {
	  m_runId   = 0;
	  m_eventId = 0;

	  memset(m_calXtalAdcPed,0,sizeof(m_calXtalAdcPed));
	  memset(m_calXtalAdcPed,0,sizeof(m_calXtalAdcPed));
  }

  int m_runId;
  int m_eventId;

  /// ped subtracted adcs
  Float_t m_calXtalAdcPed[16][8][12][2];
        
  /// Cal Xtal Face signal in scintillated MeV units.
  Float_t m_calXtalFaceSignal[16][8][12][2]; 
};


