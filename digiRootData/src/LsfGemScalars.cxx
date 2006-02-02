#include "digiRootData/LsfGemScalars.h"
#include <iostream>

ClassImp(GemScalars) ;

void GemScalars::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "GemScalars FIXME" << endl; 
}

void GemScalars::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t GemScalars::CompareInRange( const GemScalars& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
