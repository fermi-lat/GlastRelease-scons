
#ifndef AcdCalibObj_h
#define AcdCalibObj_h

// local includes
#include "CalibData/Acd/AcdCalibEnum.h"

// stl includes
#include <map>
#include <vector>
#include <string>
#include <iostream>

// forward declares


namespace CalibData {


  class AcdCalibDescription {    
  public:
    
    /// get the description for a particular calibration
    static const AcdCalibDescription* getDesc(AcdCalibData::CALTYPE calType, int version=-1) {
      const std::vector< const AcdCalibDescription* >& descs =  s_descs[calType];
      return version < 0 ? descs.back() : descs[version];
    }
  protected:
    
    /// add the description for a particular calibration type
    /// returns the version number of this particular calibration description
    static int addDesc(AcdCalibData::CALTYPE calType, const AcdCalibDescription* desc);
    
  private:
    
    /// all of the descriptions
    static std::vector< std::vector<const AcdCalibDescription*> > s_descs;
    
  public:
    
    /// Make a description
    AcdCalibDescription(AcdCalibData::CALTYPE calibType, std::string calibTypeName);
    virtual ~AcdCalibDescription() {;}
    
    inline const std::string& calibTypeName() const {
      return m_calibTypeName;
    }
    inline AcdCalibData::CALTYPE calibType() const {
      return m_calibType;
    } 
    inline int getVersion() const {
      return m_version;
    }
    inline int nVar() const {
    return m_varNames.size();
    }
    int getIndex(const std::string& name) const;
    const std::string& getVarName(int i) const;
    
  protected:
    
    void addVarName(const std::string& name);
    
  private:
    
    const AcdCalibData::CALTYPE m_calibType;
    const std::string m_calibTypeName;
    std::vector<std::string> m_varNames;
    int m_version;
    
  };
  
  
  class AcdCalibObj {
    
  public:
    
    enum STATUS { NOFIT = -1,
		  OK = 0,	 
		  // 1-3 are reserved for minuit
		  MINUIT_FAILED = 4,
		  PREFIT_FAILED = 5,
		  USED_FALLBACK_1 = 6,
		  USED_FALLBACK_2 = 7 };
    
  public:
    
    AcdCalibObj(STATUS status, const std::vector<float>& val, const AcdCalibDescription& desc);
    AcdCalibObj(STATUS status, const AcdCalibDescription& desc);
   
    virtual ~AcdCalibObj(){;}
    
    void printTxtLine(std::ostream&  os, const AcdCalibDescription& desc) const;
    bool readTxt(std::istream& is, const AcdCalibDescription& desc);
    
    inline float& operator[](int i) { return m_vals[i]; }
    inline const float& operator[](int i) const { return m_vals[i]; }
        
    inline void update(AcdCalibObj* other) {
      if ( this == other ) return;
      m_status = other->getStatus();
      m_vals.resize( other->size() );
      for ( unsigned i(0); i < m_vals.size(); i++ ) m_vals[i] = (*other)[i];
    }
    
    inline unsigned int size() const {
      return m_vals.size();
    }
    
    inline void setStatus(STATUS stat) { m_status = stat; };
    inline STATUS getStatus() const { return m_status; }
    
    void setVals(const std::vector<float>& vals, STATUS status);
    void setVals(float v1, float v2, STATUS stat);
    void setVals(float v1, float v2, float v3, STATUS stat);
    void setVals(float v1, float v2, float v3, float v4, STATUS stat);
    
  protected:
    
    AcdCalibObj();
    bool setStatusInt(int stat);
    
    inline void sizeVals(int n) {
      m_vals.resize(n);
    }
    
  private:
    
    STATUS m_status;
    
    std::vector<float> m_vals;
  };
}


#endif

