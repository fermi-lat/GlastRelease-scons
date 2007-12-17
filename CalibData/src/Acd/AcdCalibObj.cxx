#define AcdCalibResult_cxx

#include "CalibData/Acd/AcdCalibObj.h"
#include <cassert>
namespace CalibData {

    class AcdCalibDescriptionData
    {
    public:
        static AcdCalibDescriptionData* instance()
        {
            if (m_instance == 0)
            {
                m_instance = new AcdCalibDescriptionData();
            }
            return m_instance;
        }
        std::vector< std::vector<const AcdCalibDescription*> > s_descs;
    private:
        AcdCalibDescriptionData() {s_descs.resize(AcdCalibData::NDESC);}

        static AcdCalibDescriptionData* m_instance;
    };

    AcdCalibDescriptionData* AcdCalibDescriptionData::m_instance = 0;

  //std::vector< std::vector<const AcdCalibDescription*> > AcdCalibDescription::s_descs(AcdCalibData::NDESC);

  /// get the description for a particular calibration
  const AcdCalibDescription* AcdCalibDescription::getDesc(AcdCalibData::CALTYPE calType, int version) {
//      const std::vector< const AcdCalibDescription* >& descs =  s_descs[calType];
      const std::vector< const AcdCalibDescription* >& descs = AcdCalibDescriptionData::instance()->s_descs[calType];
      return version < 0 ? descs.back() : descs[version];
  }
  
  int AcdCalibDescription::addDesc(AcdCalibData::CALTYPE calType, const AcdCalibDescription* desc) {
    std::vector< std::vector<const AcdCalibDescription*> >& s_descs = AcdCalibDescriptionData::instance()->s_descs;
    if ( s_descs.size() == 0 ) {
      s_descs.resize(AcdCalibData::NDESC);
    }
    s_descs[calType].push_back(desc);
    return s_descs[calType].size();
  }
  
  AcdCalibDescription::AcdCalibDescription(AcdCalibData::CALTYPE calibType, std::string calibTypeName)
    :m_calibType(calibType),
     m_calibTypeName(calibTypeName){
    m_version = addDesc(calibType,this);
  }
  
  int AcdCalibDescription::getIndex(const std::string& name) const {
    for ( unsigned i(0); i < m_varNames.size(); i++ ) {
      if ( m_varNames[i] == name ) return i;
    }
    return -1;
  }
  
  const std::string& AcdCalibDescription::getVarName(int i) const {
    static const std::string noneName("NONE");
    if ( i < 0 || i >= (int)m_varNames.size() ) return noneName;
    return m_varNames[i];
  }
  
  void AcdCalibDescription::addVarName(const std::string& name) {
    m_varNames.push_back(name);
  }

  
  AcdCalibObj::AcdCalibObj(STATUS status, const std::vector<float>& vals, const AcdCalibDescription& desc) 
    :m_status(status){
    assert(vals.size() <= desc.nVar());
    sizeVals(desc.nVar());   
    setVals(vals,status);
  }

  AcdCalibObj::AcdCalibObj(STATUS status, const AcdCalibDescription& desc) 
    :m_status(status){
    sizeVals(desc.nVar());
  }  
    
  void AcdCalibObj::setVals(float v1, float v2, STATUS stat) {
    m_vals[0] = v1;
    m_vals[1] = v2;
    setStatus(stat);
  }
  
  void AcdCalibObj::setVals(float v1, float v2, float v3, STATUS stat){
    m_vals[0] = v1;
    m_vals[1] = v2;
    m_vals[2] = v3;
    setStatus(stat);
  }
  
  void AcdCalibObj::setVals(float v1, float v2, float v3, float v4, STATUS stat){
      sizeVals(4);
    m_vals[0] = v1;
    m_vals[1] = v2;
    m_vals[2] = v3;
    m_vals[3] = v4;
    setStatus(stat);
  }
  
  void AcdCalibObj::setVals(const std::vector<float>& vals, STATUS status) {
    if ( vals.size() == m_vals.size() ) {
      for ( unsigned i(0); i < vals.size(); i++ ) {
	m_vals[i] = vals[i];
      }
    } else {    
      return;
    }
    m_status = status;
  }  
  
  void AcdCalibObj::printTxtLine(std::ostream&  os, const AcdCalibDescription& desc ) const {
    for ( int i(0); i < desc.nVar() ; i++ ) {
      os << (float)m_vals[i] << ' ';
    }
    os << getStatus();
  }
  
  bool AcdCalibObj::readTxt(std::istream& is, const AcdCalibDescription& desc ) {
    for ( int i(0); i < desc.nVar() ; i++ ) {
      float aVal;
      is >> aVal;
      operator[](i) = aVal;
    }
    int stat;
    is >> stat;
    setStatusInt(stat);
    return true;
  }
  
  bool AcdCalibObj::setStatusInt(int stat) {
    switch (stat) {
    case NOFIT: setStatus(NOFIT); break;
    case OK: setStatus(OK); break;
    case MINUIT_FAILED: setStatus(MINUIT_FAILED); break;
    case PREFIT_FAILED: setStatus(PREFIT_FAILED); break;
    case USED_FALLBACK_1: setStatus(USED_FALLBACK_1); break;
    case USED_FALLBACK_2: setStatus(USED_FALLBACK_2); break;
    default:
      return false;
    }
    return true;
  }
}
