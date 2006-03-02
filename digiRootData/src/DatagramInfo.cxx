#include "digiRootData/DatagramInfo.h"

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

std::ostream& DatagramInfo::fillStream( std::ostream& s ) const {
     s << " open:     nmodes = " << modeChanges() << ", ndgms = "
        << datagrams() << "\n"
        << " open:     action = (" << openAction() << ")\n"
        << " open:     reason = (" << openReason() << ")\n"
        << " open:     crate = (" << crate() << ")\n"
        << " open:     mode = (" << mode() << ")\n"
        << " close:    action = (" << closeAction() << ")\n"
        << " close:    reason = (" << closeReason() << ")\n";
    return s;
}
