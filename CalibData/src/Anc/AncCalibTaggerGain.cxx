// $Header$

#include "CalibData/Anc/AncCalibTaggerGain.h"
#include "AncFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  AncCalibTaggerGain::AncCalibTaggerGain(unsigned nMod, unsigned nLay,
                                       unsigned nChan) :
    AncCalibBase(nMod, nLay, nChan) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    AncTaggerGain* pGains = new AncTaggerGain[size];
    for (ix = 0; ix < size; ix++) {
      m_chans[ix] = pGains; 
      ++pGains;
    }
  }

  AncCalibTaggerGain::~AncCalibTaggerGain() {
    AncTaggerGain* pGains = dynamic_cast<AncTaggerGain* >(m_chans[0]);
    if (pGains) delete [] pGains;
  }

  const CLID& AncCalibTaggerGain::classID()  {
    return CLID_Calib_ANC_TaggerGain;
  }

  bool AncCalibTaggerGain::putChan(unsigned mod, unsigned lay, unsigned chan, 
                               RangeBase* data) {
    if (!dynamic_cast<AncTaggerGain* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return AncCalibBase::putChan(mod, lay, chan, data);
  }

}
