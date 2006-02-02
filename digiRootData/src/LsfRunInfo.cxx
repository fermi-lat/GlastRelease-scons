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
