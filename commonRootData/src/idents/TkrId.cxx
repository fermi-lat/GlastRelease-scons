#include "commonRootData/idents/TkrId.h"

ClassImp(commonRootData::TkrId)

namespace commonRootData {

  void TkrId::init(unsigned towerX, unsigned towerY, unsigned tray, bool top, 
                   int view, int ladder, int wafer) {
    m_packedId = (towerX << SHIFTTowerX) + (towerY << SHIFTTowerY) +
      (tray << SHIFTTray);
    if (top) m_packedId |= (1 << SHIFTBotTop);
    m_packedId |=   (VALIDTowerY + VALIDTowerX + 
                     VALIDTray + VALIDBotTop);
    if ((view == eMeasureX) || (view == eMeasureY) ) {
      m_packedId |= ( (view << SHIFTMeas) | VALIDMeas);
    }
    if (ladder != eNoValue) {
      m_packedId |= ( ( ladder << SHIFTLadder) | VALIDLadder);
    }
    if (wafer != eNoValue) {
      m_packedId |= ( ( wafer << SHIFTWafer) | VALIDWafer);
    }
  }
        


}
