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
  bool done = m_cache == m_pocas.end() ? true : false;
  while ( ! done ) {
    if ( m_cache->arclength() < stop ) { 
      pocas.push_back(*m_cache);
      m_cache++;
      done = (m_cache == m_pocas.end() ? true : false);
    } else { 
      done = true;
    }
  }  
  return pocas.size();
}

void AcdPocaSorter::resetArcCache() {
  m_cache = m_pocas.begin();  
}


bool AcdPocaSorter::isDone() {
  return m_cache == m_pocas.end();  
}

double AcdPocaSorter::finalArc() {
  if ( m_pocas.size() == 0 ) return 0.;
  return m_pocas.back().arclength();
}



void AcdPocaSorter::runOut() {
  std::vector<AcdPocaSorter::AcdPocaHolder> left;
  unsigned n = getPocasToArclength( 1e6,  left );
  for ( unsigned i(0); i < n; i++ ) {
    std::cout << "Left " << i << ' ' << left[i].arclength() << ' ' << left[i].type() << std::endl;
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

