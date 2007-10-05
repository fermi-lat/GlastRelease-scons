#include "digiRootData/LsfKeys.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"

ClassImp(LsfKeys) 

ClassImp(LpaKeys) 
ClassImp(LciKeys) 

void LsfKeys::Print(Option_t* /*option*/) const {
  using namespace std;
  cout << "LATC_master = 0x " << setw(8) << setfill('0') << hex 
       << m_LATC_master << endl;
  cout << "LATC_ignore = 0x" << setw(8) << setfill('0') << hex
       << m_LATC_ignore << endl;
}

void LsfKeys::Fake( Int_t ievent, Float_t /* randNum */ ) {
  m_LATC_master = ievent;
  m_LATC_ignore = ievent*2;
}

Bool_t LsfKeys::CompareInRange( const LsfKeys& ref, 
			        const std::string & name ) const {
  /// FIXME
  bool result = true ;
 result = rootdatautil::CompareInRange(LATC_master(),ref.LATC_master(),"LATC_master") && result;
 result = rootdatautil::CompareInRange(LATC_ignore(),ref.LATC_ignore(),"LATC_ignore") && result;

  return result;
}

void LpaKeys::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "LpaKeys: " << endl;
  cout << " Master = 0x" << std::hex << std::setfill('0')
       << m_LATC_master << "\n"
       << " Ignore = 0x" << m_LATC_ignore << std::endl;
  std::vector<UInt_t>::const_iterator it;
  unsigned int i = 0;
  for(it = m_CDM_keys.begin(); it != m_CDM_keys.end(); it++) {
    cout << std::dec << "CDM key " << i << " " << m_CDM_keys[i] << endl;
    ++i;
  }
}

void LpaKeys::Fake( Int_t ievent, Float_t randNum ) {
  LsfKeys::Fake(ievent, randNum);
  std::vector<unsigned int> vec;
  unsigned int i;
  for (i=0; i<10; i++) 
    vec.push_back(ievent*i);
  setCDM_keys(vec);
}

Bool_t LpaKeys::CompareInRange( const LpaKeys& ref, 
			        const std::string & name ) const {
  bool result = true ;
  result = rootdatautil::CompareInRange(CDM_keys(),ref.CDM_keys(),"CDM_keys") && result;
 result = rootdatautil::CompareInRange(LATC_master(),ref.LATC_master(),"LATC_master") && result;
 result = rootdatautil::CompareInRange(LATC_ignore(),ref.LATC_ignore(),"LATC_ignore") && result;
  return result;
}


void LciKeys::Print(Option_t* /*option*/) const {
  /// FIXME
  using namespace std;
  cout << "LciKeys: " <<  endl;
  cout << "Master = 0x" << std::hex << std::setfill('0')
       << m_LATC_master << endl
       << "Ignore: " << m_LATC_ignore << endl
       << "LCI_Script: " << m_LCI_script << endl;
}

void LciKeys::Fake( Int_t ievent, Float_t randNum ) {
  LsfKeys::Fake(ievent, randNum);
  setLCI_script(ievent);
}

Bool_t LciKeys::CompareInRange( const LciKeys& ref, 
        const std::string & name ) const {
  bool result = true ;
  result = rootdatautil::CompareInRange(LCI_script(),ref.LCI_script(),"LCI_script") && result;
 result = rootdatautil::CompareInRange(LATC_master(),ref.LATC_master(),"LATC_master") && result;
 result = rootdatautil::CompareInRange(LATC_ignore(),ref.LATC_ignore(),"LATC_ignore") && result;
  return result;
}

