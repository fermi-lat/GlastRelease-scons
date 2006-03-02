#include "digiRootData/LsfTime.h"

ClassImp(LsfTime) ;

void LsfTime::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "Time FIXME" << endl; 
}

std::ostream& LsfTime::fillStream( std::ostream& s ) const {
      s << " current:  secs = " << m_current.timeSecs() << "\n"
        << "   tics  = 0x" << std::hex << std::uppercase << std::setfill('0')
        << std::setw(8)
        << m_current.timeHack().ticks() << " ("
        << std::dec << m_current.timeHack().ticks() << ")" << std::endl
        << "   hacks = 0x" << std::hex << m_current.timeHack().hacks() << " ("
        << std::dec << m_current.timeHack().hacks() << ")" << std::endl
        << " previous:  secs = " << m_previous.timeSecs() << "\n"
        << "   tics  = 0x" << std::hex << m_previous.timeHack().ticks() << " ("
        << std::dec << m_previous.timeHack().ticks() << ")" << std::endl
        << "   hacks = 0x" << std::hex << m_previous.timeHack().hacks() << " ("
        << std::dec << m_previous.timeHack().hacks() << ")" << std::endl;
    return s;
}

void LsfTime::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LsfTime::CompareInRange( const LsfTime& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
