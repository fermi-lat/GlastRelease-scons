#include "digiRootData/GemScalers.h"

ClassImp(GemScalers) ;

void GemScalers::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "GemScalers FIXME" << endl; 
}

std::ostream& GemScalers::fillStream( std::ostream& s ) const {
      s << " scalers:  elapsed   = 0x" << std::hex << std::nouppercase
        << std::setfill('0') << std::setw(8)
        << elapsed() << " = " << std::dec << elapsed() << "\n"
        << " scalers:  livetime  = 0x" << std::hex << livetime() << " = "
        << std::dec << livetime() << "\n"
        << " scalers:  prescaled = 0x" << std::hex << prescaled() << " = "
        << std::dec << prescaled() << "\n"
        << " scalers:  discarded = 0x" << std::hex << discarded() << " = "
        << std::dec << discarded() << "\n"
        << " scalers:  deadzone  = 0x" << std::hex << deadzone() << " = "
        << std::dec << deadzone() << "\n";
    return s;
}

void GemScalers::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t GemScalers::CompareInRange( const GemScalers& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
