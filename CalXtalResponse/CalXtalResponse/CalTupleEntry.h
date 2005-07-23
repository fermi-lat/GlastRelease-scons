// LOCAL INCLUDES

// GLAST INCLUDES

// EXTLIB INCLUDES
#include "Rtypes.h"

// STD INCLUDES
#include <cstring>

/// Single entry in CalTuple
struct CalTupleEntry {
  void Clear() {memset(this,0,sizeof(CalTupleEntry));}

  int m_runId;
  int m_eventId;

  /// ped subtracted adcs
  Float_t m_calXtalAdcPed[16][8][12][2];
        
  /// Cal Xtal Face signal in scintillated MeV units.
  Float_t m_calXtalFaceSignal[16][8][12][2];

  static const char *m_tupleDesc;
};


