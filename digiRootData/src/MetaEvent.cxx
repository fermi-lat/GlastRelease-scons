#include "digiRootData/MetaEvent.h"

ClassImp(MetaEvent) ;

void MetaEvent::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << m_run << m_datagram << endl;
  cout << "=========================================="
            << std::endl << std::endl << "Event "
            << m_scalers.sequence() << " context:" << std::endl;

  cout << m_time << m_scalers << m_run << m_datagram << std::endl
       << "Event " << m_scalers.sequence() << " info:" << std::endl
       << "---------------" << std::endl
       << "LSE_Info:  timeTicks = 0x" << std::uppercase << std::hex
       << std::setfill('0') << std::setw(8)
       << m_time.timeTicks() << " (" << std::dec
       << m_time.timeTicks()
       << ")\nLSE_Info:  " << m_time.timeHack() << std::endl;
      if (m_config) {
          const LpaConfiguration* lpa = m_config->castToLpaConfig();
          if (lpa) {
              cout << "LPA_Info: softwareKey = 0x" << std::hex
                   << lpa->softwareKey() << std::endl
                   << "LPA_Info: hardwareKey = 0x" << lpa->hardwareKey()
                   << std::dec << std::endl;
          }
      }
}

void MetaEvent::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t MetaEvent::CompareInRange( const MetaEvent& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
