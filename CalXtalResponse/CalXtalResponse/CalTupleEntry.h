// LOCAL INCLUDES

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES
#include <cstring>

/// Single entry in CalTuple
struct CalTupleEntry {
  void Clear() {memset(this,0,sizeof(CalTupleEntry));}

  int m_runId;
  int m_eventId;

  /// ped subtracted adcs
  float m_calXtalAdcPed[16][8][12][2];
        
  /// Cal Xtal Face signal in scintillated MeV units.
  float m_calXtalFaceSignal[16][8][12][2];
};


