#include "digiRootData/Ccsds.h"

ClassImp(Ccsds) ;

void Ccsds::Print(Option_t* /*option*/) const {
  using namespace std;
  cout << "CCSDS " << std::endl << "SCID: " <<  m_scid 
       << " APID: " << m_apid << " UTC: " << m_utc << std::endl;
}

void Ccsds::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t Ccsds::CompareInRange( const Ccsds& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
