
#include "commonRootData/idents/TkrId.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

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
        
void TkrId::Fake( Int_t /* ievent */, UInt_t rank, Float_t /* randNum */ ) {

    UInt_t towerX = rank % 4 ;
    UInt_t towerY = ( rank / 4 ) % 4 ;
    UInt_t tray   = rank % 15 ;
    init(towerX,towerY,tray,true,eMeasureX ) ;

}

#define COMPARE_IN_RANGE(att,text) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),text)

Bool_t TkrId::CompareInRange( const TkrId & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(View,"Axes") && result ;
    result = COMPARE_IN_RANGE(TowerX,"TowerX") && result ;
    result = COMPARE_IN_RANGE(TowerY,"TowerY") && result ;
    result = isEqual(ref) && result ; // this single one would be enough !

    if (!result) {
        std::cout<<"Comparison ERROR for "<<name<<std::endl ;
    }
    return result ;

}

}
