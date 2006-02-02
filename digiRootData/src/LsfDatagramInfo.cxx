#include "digiRootData/LsfDatagramInfo.h"
#include <iostream>

ClassImp(DatagramInfo) ;

void DatagramInfo::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "DatagramInfo FIXME" << endl; 
}

void DatagramInfo::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t DatagramInfo::CompareInRange( const DatagramInfo& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
