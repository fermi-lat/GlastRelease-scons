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

// Constructor from VolumeIdentifier.  Assuming a particular form for
// fields within VolumeIdentifier (correct as of June 14, 2004):
//    field 0 is fLATObjects; check for value = 0 (towers)
//    field 1 is y-value for tower 
//    field 2 is x-value for tower
//    field 3 is fTowerObjects; check for value = 1 (Tracker)
//    field 4 is tray number
//    field 5 is orientation (measures X or Y)
//    field 6 is top or bottom (Si for tray)
//    field 7 is ladder number
//    field 8 is wafer number
TkrId::TkrId(const VolumeIdentifier& vId) : m_packedId(0), m_validFields(0) {
  constructorGuts(vId);
}

// the inserter; expect at most diagnostic use
void TkrId::write(std::ostream &stream) const
{
  stream << "(" << std::hex << m_packedId << " ,";
  stream << m_validFields << ")";
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
  if ((vId[fTowerY] > 4) || (vId[fTowerX] > 4)) 
    throw std::invalid_argument("VolumeIdentifier");
  m_packedId = (vId[fTowerY] << SHIFTTowerY) + (vId[fTowerX] << SHIFTTowerX);
  m_validFields = SHMASKTowerY + SHMASKTowerX;
  
  if (vIdSize > fTray) {
    m_packedId |= (vId[fTray] << SHIFTTray);
    m_validFields |= SHMASKTray;
  }
  
  if (vIdSize > fMeas) {
    m_packedId |= (vId[fMeas] << SHIFTMeas);
    m_validFields |= SHMASKMeas;
  }
  
  if (vIdSize > fBotTop) {
    m_packedId |= (vId[fBotTop] << SHIFTBotTop);
    m_validFields |= SHMASKBotTop;
    
    if (vIdSize > fLadder) {
      m_packedId |= (vId[fLadder] << SHIFTLadder);
      m_validFields |= SHMASKLadder;
    }
    
    if (vIdSize > fWafer) {
      m_packedId |= (vId[fWafer] << SHIFTWafer);
      m_validFields |= SHMASKWafer;
    }    
  }
}





