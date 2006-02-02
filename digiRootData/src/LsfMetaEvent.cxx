#include "digiRootData/LsfMetaEvent.h"
#include <iostream>

ClassImp(MetaEvent) ;

void MetaEvent::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "MetaEvent FIXME" << endl; 
}

void MetaEvent::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t MetaEvent::CompareInRange( const MetaEvent& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
