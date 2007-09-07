#include "./AcdPocaSorter.h"
#include <cassert>
#include <algorithm>


double AcdPocaSorter::AcdPocaHolder::arclength() const {
  if ( m_poca == 0 ) return 0.;
  switch( m_type ) {
  case PlanePoca: return m_poca->m_arcLengthPlane;
  case RayPoca: return m_poca->m_arcLength;
  }
  return 0.;
}

bool AcdPocaSorter::AcdPocaHolder::operator==(const AcdPocaHolder& other) const {
  if ( this == &other) return true;
  if ( type() != other.type() ) return false;
  if ( arclength() != other.arclength() ) return false;
  if ( m_poca->m_id.id() != other.poca()->m_id.id() ) return false;
  assert(0);
}

bool AcdPocaSorter::AcdPocaHolder::operator<(const AcdPocaHolder& other) const{
  if ( arclength() < other.arclength() ) return true;
  if ( arclength() > other.arclength() ) return false;
  if ( type() < other.type() ) return true;
  if ( type() > other.type() ) return false;
  if ( m_poca->m_id.id() < other.poca()->m_id.id() ) return true;
  if ( m_poca->m_id.id() > other.poca()->m_id.id() ) return false;
  return false;  
}

bool AcdPocaSorter::AcdPocaHolder::operator>(const AcdPocaHolder& other) const{
  if ( arclength() > other.arclength() ) return true;
  if ( arclength() < other.arclength() ) return false;
  if ( type() > other.type() ) return true;
  if ( type() < other.type() ) return false;
  if ( m_poca->m_id.id() > other.poca()->m_id.id() ) return true;
  if ( m_poca->m_id.id() < other.poca()->m_id.id() ) return false;
  return false;   
}

bool AcdPocaSorter::AcdPocaHolder::operator<=(const AcdPocaHolder& other) const{
  if ( this == &other ) return true;
  if ( arclength() < other.arclength() ) return true;
  if ( arclength() > other.arclength() ) return false;
  if ( type() < other.type() ) return true;
  if ( type() > other.type() ) return false;
  if ( m_poca->m_id.id() < other.poca()->m_id.id() ) return true;
  if ( m_poca->m_id.id() > other.poca()->m_id.id() ) return false;
  assert(0);
}

bool AcdPocaSorter::AcdPocaHolder::operator>=(const AcdPocaHolder& other) const{
  if ( this == &other ) return true;
  if ( arclength() > other.arclength() ) return true;
  if ( arclength() < other.arclength() ) return false;
  if ( type() > other.type() ) return true;
  if ( type() < other.type() ) return false;
  if ( m_poca->m_id.id() > other.poca()->m_id.id() ) return true;
  if ( m_poca->m_id.id() < other.poca()->m_id.id() ) return false;
  assert(0);
}

AcdPocaSorter::AcdPocaSorter(TrackDirection dir, const AcdRecon::PocaDataPtrMap& theMap)
  :m_dir(dir){
  for (  AcdRecon::PocaDataPtrMap::const_iterator itr = theMap.begin(); itr != theMap.end(); itr++ ) {
    addPoca(itr->second);
  }
  sort();
}



unsigned AcdPocaSorter::getPocasToArclength(const double& stop, std::vector<AcdPocaSorter::AcdPocaHolder>& pocas) {
  pocas.clear();
  bool done(false);
  switch ( m_dir ) {
  case Upward: done = (m_cache == m_pocas.end() ? true : false);  break;    
  case Downward: done = (m_rcache == m_pocas.rend() ? true : false);  break;
  }  

  while ( ! done ) {
    switch ( m_dir ) {
    case Upward: 
      if ( m_cache->arclength() < stop ) { 
	pocas.push_back(*m_cache);
	m_cache++;
	done = (m_cache == m_pocas.end() ? true : false);
      } else { 
	done = true;
      }
      break;
    case Downward: 
      if ( m_rcache->arclength() > stop ) {
	pocas.push_back(*m_rcache);
	m_rcache++;
	done = (m_rcache == m_pocas.rend() ? true : false);
      } else {
	done = true;
      }
      break;
    default:
      break;
    }
  }  
  return pocas.size();
}

void AcdPocaSorter::resetArcCache() {
  switch ( m_dir ) {
  case Upward:  m_cache = m_pocas.begin();  break;
  case Downward:  m_rcache = m_pocas.rbegin();  break;
  }        
}


void AcdPocaSorter::sort() {
    std::sort(m_pocas.begin(),m_pocas.end());
  resetArcCache();
}

void AcdPocaSorter::addPoca(const AcdRecon::PocaData* pocaData) {
  AcdRecon::PocaData* ncP = const_cast<AcdRecon::PocaData*>(pocaData);
  m_pocas.push_back( AcdPocaSorter::AcdPocaHolder(AcdPocaSorter::PlanePoca,ncP));
  m_pocas.push_back( AcdPocaSorter::AcdPocaHolder(AcdPocaSorter::RayPoca,ncP));  
}
