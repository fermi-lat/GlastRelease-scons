#include "digiRootData/Configuration.h"

ClassImp(Configuration) 

ClassImp(LpaConfiguration) 
ClassImp(LciConfiguration) 
ClassImp(LciAcdConfiguration) 
ClassImp(LciCalConfiguration) 
ClassImp(LciTkrConfiguration) 
ClassImp(LciAcdConfiguration::AcdTrigger)
ClassImp(LciCalConfiguration::CalTrigger)

void LpaConfiguration::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << " softwareKey = 0x" << std::hex << std::setfill('0')
       << m_softwareKey << "\n"
       << " hardwareKey = 0x" << m_hardwareKey << std::endl;
}

void LpaConfiguration::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LpaConfiguration::CompareInRange( const LpaConfiguration& /* other */, 
					 const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}


void LciConfiguration::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "softwareKey = 0x" << std::hex << std::setfill('0')
       << m_softwareKey << endl << std::dec << "writeCfg: " << m_writeCfg 
       << endl
       << "readCfg: " << m_readCfg << endl << "period: " << m_period << endl
       << "flags: " << m_flags << endl;
}

void LciConfiguration::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LciConfiguration::CompareInRange( const LciConfiguration& /* other */, 
        const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

void Channel::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t Channel::CompareInRange( const Channel& /* other */, 
       const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

void LciAcdConfiguration::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LciAcdConfiguration::CompareInRange( const LciAcdConfiguration& /* other */, 
   const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

void LciAcdConfiguration::AcdTrigger::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LciAcdConfiguration::AcdTrigger::CompareInRange( 
             const LciAcdConfiguration::AcdTrigger& /* other */, 
             const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

void LciCalConfiguration::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LciCalConfiguration::CompareInRange( 
      const LciCalConfiguration& /* other */, 
      const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}

void LciCalConfiguration::CalTrigger::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LciCalConfiguration::CalTrigger::CompareInRange( 
             const LciCalConfiguration::CalTrigger& /* other */, 
             const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}


void LciTkrConfiguration::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  /// FIXME
}

Bool_t LciTkrConfiguration::CompareInRange( 
   const LciTkrConfiguration& /* other */, 
   const std::string & /* name = "" */ ) const {
  /// FIXME
  return kTRUE;
}
