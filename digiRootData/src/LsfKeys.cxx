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
       << " Ignore = 0x" << m_LATC_ignore 
       << " SBS: " << m_sbs << " LPA_DB: " << m_lpa_db << std::endl;
}

void LpaKeys::Fake( Int_t ievent, Float_t randNum ) {
  LsfKeys::Fake(ievent, randNum);
  m_sbs = ievent * 3;
  m_lpa_db = ievent * 4;
}

Bool_t LpaKeys::CompareInRange( const LpaKeys& ref, 
			        const std::string & name ) const {
  bool result = true ;
  result = rootdatautil::CompareInRange(sbs(),ref.sbs(),"SBS") && result;
  result = rootdatautil::CompareInRange(lpa_db(),ref.lpa_db(),"LPA_DB") && result;
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

