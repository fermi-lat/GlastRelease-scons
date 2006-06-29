// $Header$

#include "CalibData/Anc/AncCalibQdcPed.h"
#include "AncFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AncCalibQdcPed::AncCalibQdcPed(unsigned nMod, unsigned nChan) :
    AncCalibBase(nMod, 1, nChan) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AncQdcPed* pPeds = new AncQdcPed[size];
    for (ix = 0; ix < size; ix++) {
      m_chans[ix] = pPeds; 
      ++pPeds;
    }
  }

  AncCalibQdcPed::~AncCalibQdcPed() {
    AncQdcPed* pPeds = dynamic_cast<AncQdcPed* >(m_chans[0]);
    if (pPeds) delete [] pPeds;
  }

  const CLID& AncCalibQdcPed::classID()  {return CLID_Calib_ANC_QdcPed;}

  bool AncCalibQdcPed::putChan(unsigned mod, unsigned lay, unsigned chan, 
                               RangeBase* data) {
    if (!dynamic_cast<AncQdcPed* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AncCalibBase::putChan(mod, lay, chan, data);
  }

}
