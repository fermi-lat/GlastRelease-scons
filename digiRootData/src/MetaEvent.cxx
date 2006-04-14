#include "digiRootData/MetaEvent.h"

ClassImp(MetaEvent) ;

void MetaEvent::Print(Option_t* /*option*/) const {
  using namespace std;
  cout << "=========================================="
            << std::endl << std::endl << "Event "
            << m_scalers.sequence() << " context:" << std::endl;

       m_time.Print("");
       m_scalers.Print("");
       m_run.Print("");
       m_datagram.Print("");
       cout << "Event " << m_scalers.sequence() << " info:" << std::endl
       << "---------------" << std::endl
       << "LSE_Info:  timeTicks = 0x" << std::uppercase << std::hex
       << std::setfill('0') << std::setw(8)
       << m_time.timeTicks() << " (" << std::dec
       << m_time.timeTicks()
       << ")\nLSE_Info:  ";
       m_time.timeHack().Print("");
       cout << std::endl;
/*
      if (m_config) {
          const LpaConfiguration* lpa = m_config->castToLpaConfig();
          if (lpa) {
              cout << "LPA_Info: softwareKey = 0x" << std::hex
                   << lpa->softwareKey() << std::endl
                   << "LPA_Info: hardwareKey = 0x" << lpa->hardwareKey()
                   << std::dec << std::endl;
          }
      }
*/
      if (lpaConfiguration())
         lpaConfiguration()->Print("");
      if (lciAcdConfiguration())
         lciAcdConfiguration()->Print("");
      if (lciCalConfiguration())
         lciCalConfiguration()->Print("");
      if (lciTkrConfiguration())
         lciTkrConfiguration()->Print("");
}

void MetaEvent::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t MetaEvent::CompareInRange( const MetaEvent& /* other */, const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
