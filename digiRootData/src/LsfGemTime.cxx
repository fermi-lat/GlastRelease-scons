#include "digiRootData/LsfGemTime.h"
#include <iostream>

ClassImp(GemTime) ;

void GemTime::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "GemTime FIXME" << endl; 
}

void GemTime::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t GemTime::CompareInRange( const GemTime& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

