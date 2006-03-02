#include "digiRootData/GemTime.h"

ClassImp(GemTime) ;

void GemTime::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "GemTime FIXME" << endl; 
}

std::ostream& GemTime::fillStream( std::ostream& s ) const {
    return s;
}

void GemTime::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t GemTime::CompareInRange( const GemTime& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

