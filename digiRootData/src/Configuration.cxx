#include "digiRootData/Configuration.h"
#include <iostream>

ClassImp(Configuration) ;

ClassImp(LpaConfiguration) ;

void LpaConfiguration::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "Configuration FIXME" << endl; 
}

void LpaConfiguration::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LpaConfiguration::CompareInRange( const LpaConfiguration& /* other */, 
					 const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}


ClassImp(LciConfiguration) ;

void LciConfiguration::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "Configuration FIXME" << endl; 
}

void LciConfiguration::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LciConfiguration::CompareInRange( const LciConfiguration& /* other */, 
					 const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
