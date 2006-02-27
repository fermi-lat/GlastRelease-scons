#include "digiRootData/LsfTime.h"
#include <iostream>

ClassImp(LsfTime) ;

void LsfTime::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "Time FIXME" << endl; 
}

void LsfTime::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LsfTime::CompareInRange( const LsfTime& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
