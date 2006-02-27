#include "digiRootData/GemScalers.h"
#include <iostream>

ClassImp(GemScalers) ;

void GemScalers::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "GemScalers FIXME" << endl; 
}

void GemScalers::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t GemScalers::CompareInRange( const GemScalers& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
