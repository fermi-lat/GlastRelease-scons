// File and Version information
// $Header$
//
// ClassName:   TkrId
//  
// Description: ID class for Tracker components
//
//              
// Author:  J. Bogart     14 June 2004

// Include files
#include "idents/TkrId.h"
#include "idents/VolumeIdentifier.h"
#include <stdexcept>
#include <iostream>
#include <ios>

using namespace idents; 

/** Constructor from VolumeIdentifier.  Assuming a particular form for
    fields within VolumeIdentifier (correct as of June 14, 2004):
    field 0 is fLATObjects; check for value = 0 (towers)
    field 1 is y-value for tower 
    field 2 is x-value for tower
    field 3 is fTowerObjects; check for value = 1 (Tracker)
    field 4 is tray number
    field 5 is orientation (measures X or Y)
    field 6 is top or bottom (Si for tray)
    field 7 is ladder number
    field 8 is wafer number
*/
TkrId::TkrId(const VolumeIdentifier& vId) : m_packedId(0) {
  constructorGuts(vId);
}

/** 
 Simple alternate constructor when all we care about is the plane. 
 Designed for use with calibration data or, with optional view
 argument, analysis code.
*/
TkrId::TkrId(unsigned towerX, unsigned towerY, unsigned tray, bool top,
             int view) {
  m_packedId = (towerX << SHIFTTowerX) + (towerY << SHIFTTowerY) +
    (tray << SHIFTTray);
  if (top) m_packedId |= (1 << SHIFTBotTop);
  m_packedId |=   (VALIDTowerY + VALIDTowerX + 
                   VALIDTray + VALIDBotTop);
  if ((view == eMeasureX) || (view == eMeasureY) ) {
    m_packedId |= ( (view << SHIFTMeas) | VALIDMeas);
  }
}

// the inserter; expect at most diagnostic use
void TkrId::write(std::ostream &stream) const
{
  stream << "(" << std::hex << m_packedId << " )";
  stream << std::dec ;          // return to presumed normal state
}

void TkrId::constructorGuts(const VolumeIdentifier& vId)   {
  const int minSize = 4;
  const unsigned eLATTowers = 0;
  const unsigned eTowerTKR = 1;
  enum {
    fLATObjects = 0,
    fTowerY = 1,
    fTowerX = 2,
    fTowerObjects = 3,
    fTray = 4,
    fMeas = 5, 
    fBotTop = 6,
    fLadder = 7,
    fWafer = 8
  };
  int vIdSize = vId.size();
  if (vIdSize < minSize) throw std::invalid_argument("VolumeIdentifier");
  if ((vId[fLATObjects] != eLATTowers) || 
      (vId[fTowerObjects] != eTowerTKR)) {
    throw std::invalid_argument("VolumeIdentifier");
  }
  if ((vId[fTowerY] > 3) || (vId[fTowerX] > 3)) 
    throw std::invalid_argument("VolumeIdentifier");
  m_packedId = (vId[fTowerY] << SHIFTTowerY) + (vId[fTowerX] << SHIFTTowerX);
  m_packedId |= VALIDTowerY + VALIDTowerX;
  
  if (vIdSize > fTray) {
    m_packedId |= (vId[fTray] << SHIFTTray);
    m_packedId |= VALIDTray;
  }
  
  if (vIdSize > fMeas) {
    m_packedId |= (vId[fMeas] << SHIFTMeas);
    m_packedId |= VALIDMeas;
  }
  
  if (vIdSize > fBotTop) {
    m_packedId |= (vId[fBotTop] << SHIFTBotTop);
    m_packedId |= VALIDBotTop;
    
    if (vIdSize > fLadder) {
      m_packedId |= (vId[fLadder] << SHIFTLadder);
      m_packedId |= VALIDLadder;
    }
    
    if (vIdSize > fWafer) {
      m_packedId |= (vId[fWafer] << SHIFTWafer);
      m_packedId |= VALIDWafer;
    }    
  }
}





