#include "McHitsManager.h"

// Gaudi
#include "GaudiKernel/SmartDataPtr.h"


// GlastEvent for creating the McEvent stuff
#include "GlastEvent/MonteCarlo/McPositionHit.h"

McHitsManager* McHitsManager::pointer = 0;

McHitsManager* McHitsManager::getPointer()
{
  if(pointer == 0)
    pointer = new McHitsManager;
  return pointer;
}

void McHitsManager::addHit(mc::McPositionHit *hit)
{
  m_posHit->push_back(hit);
}







