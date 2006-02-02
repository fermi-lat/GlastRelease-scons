#include "digiRootData/LsfTimeHack.h"
#include <iostream>

ClassImp(TimeHack) ;

void TimeHack::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "TimeHack FIXME" << endl; 
}

void TimeHack::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t TimeHack::CompareInRange( const TimeHack& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

