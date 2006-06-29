// $Header$

#include "CalibData/Anc/AncCalibTaggerPed.h"
#include "AncFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AncCalibTaggerPed::AncCalibTaggerPed(unsigned nMod, unsigned nLay,
                                       unsigned nChan) :
    AncCalibBase(nMod, nLay, nChan) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AncTaggerPed* pPeds = new AncTaggerPed[size];
    for (ix = 0; ix < size; ix++) {
      m_chans[ix] = pPeds; 
      ++pPeds;
    }
  }

  AncCalibTaggerPed::~AncCalibTaggerPed() {
    AncTaggerPed* pPeds = dynamic_cast<AncTaggerPed* >(m_chans[0]);
    if (pPeds) delete [] pPeds;
  }

  const CLID& AncCalibTaggerPed::classID()  {return CLID_Calib_ANC_TaggerPed;}

  bool AncCalibTaggerPed::putChan(unsigned mod, unsigned lay, unsigned chan, 
                               RangeBase* data) {
    if (!dynamic_cast<AncTaggerPed* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AncCalibBase::putChan(mod, lay, chan, data);
  }

}
