// File and Version information:
// $Header$
//
//  Implementation file of AcdTkrPoca and AcdTkrPocaCol classes
//  
// Authors:
//
//    Eric Charles
//
//

#include "Event/Recon/AcdRecon/AcdPocaMap.h"

#include "GaudiKernel/MsgStream.h"


using namespace Event;

/// default constructor, builds any empty map
AcdPocaMap::AcdPocaMap(){
  ini();
}

/// adds a single poca to the map
void AcdPocaMap::add(const Event::AcdTkrHitPoca& poca) {

  int trackIdx = poca.trackIndex();
  const idents::AcdId& id = poca.getId();

  m_tileToPocaMap[id].insert(&poca);
  m_trackToPocaMap[trackIdx].insert(&poca);
}

void AcdPocaMap::writeOut(MsgStream& /* stream */ ) const

// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

}

/// get all the pocas associated with a single tile
const Event::AcdPocaSet& AcdPocaMap::getPocas(const idents::AcdId& acdId) const {
  static const AcdPocaSet nullSet;
  std::map< idents::AcdId, Event::AcdPocaSet >::const_iterator it = m_tileToPocaMap.find(acdId);
  if ( it == m_tileToPocaMap.end() ) return nullSet;
  return it->second;
}

/// get the best poca associated with a single tile (can return null)
const Event::AcdTkrHitPoca* AcdPocaMap::getBestPoca(const idents::AcdId& acdId) const {
  std::map< idents::AcdId, Event::AcdPocaSet >::const_iterator it = m_tileToPocaMap.find(acdId);
  if ( it == m_tileToPocaMap.end() ) return 0;
  Event::AcdPocaSet::const_iterator itS = it->second.begin();
  if ( itS == it->second.end() ) return 0;
  return *itS;
}

/// get all the pocas associated with a single track
const Event::AcdPocaSet& AcdPocaMap::getPocas(int trackIdx) const {
  static const AcdPocaSet nullSet;
  std::map< int, Event::AcdPocaSet >::const_iterator it = m_trackToPocaMap.find(trackIdx);
  if ( it == m_trackToPocaMap.end() ) return nullSet;
  return it->second;
}

/// get the best poca associated with a single track (can return null)
const Event::AcdTkrHitPoca* AcdPocaMap::getBestPoca(int trackIdx) const {
  std::map< int, Event::AcdPocaSet >::const_iterator it = m_trackToPocaMap.find(trackIdx);
  if ( it == m_trackToPocaMap.end() ) return 0;
  Event::AcdPocaSet::const_iterator itS = it->second.begin();
  if ( itS == it->second.end() ) return 0;
  return *itS;
}

/// get the poca between a track and a tile (can return null)
const Event::AcdTkrHitPoca* AcdPocaMap::getPoca(const idents::AcdId& acdId, int trackIdx ) const {
  const AcdPocaSet& set = getPocas(acdId);
  for ( AcdPocaSet::const_iterator it = set.begin(); it != set.end(); it++ ) {
    if ( (*it)->trackIndex() == trackIdx ) return *it;
  }
  return 0;
}

void AcdPocaMap::ini()
// Purpose: reset all data members to 0
//
{
  m_tileToPocaMap.clear();
  m_trackToPocaMap.clear();
}
