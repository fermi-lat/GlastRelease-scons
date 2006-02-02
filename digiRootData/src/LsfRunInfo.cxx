#include "digiRootData/LsfRunInfo.h"
#include <iostream>

ClassImp(RunInfo) ;

void RunInfo::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "RunInfo FIXME" << endl; 
}

void RunInfo::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t RunInfo::CompareInRange( const RunInfo& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

/// For sorting, override TObject::Compare()
Int_t RunInfo::Compare(const TObject *obj) const {
  if (this == obj) return 0;
  if (obj == 0 ) return 1;
  
  // This logic to handle class mismatches is copied from TObjString
  // Please note that this is not a good idea as it breaks 
  // strict weak ordering & can cause all sort of problems in sorted collections
  if (RunInfo::Class() != obj->IsA()) return -1;

  // cast the other as a RunInfo
  const RunInfo* other = static_cast<const RunInfo*>(obj);
  
  // compare start time first
  if ( m_startTime != other->startTime() ) return m_startTime > other->startTime() ? 1 : -1;
  // then check ground id
  if ( m_id != other->id() ) return m_id > other->id() ? 1 : -1;
  // then check platform
  if ( m_platform != other->platform() ) return m_platform > other->platform() ? 1 : -1;
  // then check origin
  if ( m_origin != other->dataOrigin() ) return m_origin > other->dataOrigin() ? 1 : -1;
  // they are equivalent
  return 0;  
}
