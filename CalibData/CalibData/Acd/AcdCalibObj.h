
#ifndef AcdCalibObj_h
#define AcdCalibObj_h

// local includes
#include "CalibData/Acd/AcdCalibEnum.h"

// stl includes
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

// forward declares

namespace CalibData {

  /** 
   * @class AcdCalibDescription
   *
   * @brief Abstract Description of an ACD calibration.
   * 
   * This class really just maps variable names to index
   * Each ACD calibration is just a few floats, this class allows use to keep track of 
   * what name goes with each of those floats.
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdCalibDescription {    
  public:
    
    /// get the description for a particular calibration
    static const AcdCalibDescription* getDesc(AcdCalibData::CALTYPE calType, 
                                              int version=-1);
  protected:
    
    /// add the description for a particular calibration type
    /// returns the version number of this particular calibration description
    static int addDesc(AcdCalibData::CALTYPE calType, 
                       const AcdCalibDescription* desc);
    
  public:
    
    /// Make a description
    AcdCalibDescription(AcdCalibData::CALTYPE calibType, 
                        std::string calibTypeName);

    /// Trivial D'ctor
    virtual ~AcdCalibDescription() {;}
    
    /// Return the name of the calibration type
    inline const std::string& calibTypeName() const {
      return m_calibTypeName;
    }
    /// Return an enum with the calibration type
    inline AcdCalibData::CALTYPE calibType() const {
      return m_calibType;
    } 
    /// Get the verions of the calibration 
    inline int getVersion() const {
      return m_version;
    }
    /// return the number of varibles in this calibration
    inline int nVar() const {
      return m_varNames.size();
    }
    /// get the index of a particular variable by name
    int getIndex(const std::string& name) const;
    /// get the name of a particular variable by index
    const std::string& getVarName(int i) const;
    
  protected:
    
    /// Add a variable to this calibration type
    void addVarName(const std::string& name);
    
  private:
    
    /// The type of calibration
    const AcdCalibData::CALTYPE m_calibType;
    /// The name of the calibration
    const std::string m_calibTypeName;
    /// The names of the variables
    std::vector<std::string> m_varNames;
    /// The version of the calibration
    int m_version;
    
  };
  
   /** 
   * @class AcdCalibObj
   *
   * @brief Base class for all calibrations of individual PMTs.
   *
   * All this class really does is associate a calibration description (AcdCalibDescription)
   * with a particular set of values.  The descrption associates names with the values by index.
   * 
   * Sub-classes provide access functions to individual value by name.
   * 
   * @author Eric Charles
   */
 
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
    
    /// Build calibration and fill values
    AcdCalibObj(STATUS status, const std::vector<float>& val, 
                const AcdCalibDescription& desc);

    /// Build calibration, but don't fill values
    AcdCalibObj(STATUS status, const AcdCalibDescription& desc);
    
    /// trivial d'tor
    virtual ~AcdCalibObj(){;}
    
    void printTxtLine(std::ostream&  os, 
                      const AcdCalibDescription& desc) const;
    bool readTxt(std::istream& is, const AcdCalibDescription& desc);

    // access a particular variable by index
    inline float& operator[](int i) { return m_vals[i]; }
    inline const float& operator[](int i) const { return m_vals[i]; }
        
    /// copy another calibration into this one
    inline void update(AcdCalibObj* other) {
      if ( this == other ) return;
      m_status = other->getStatus();
      m_vals.resize( other->size() );
      for ( unsigned i(0); i < m_vals.size(); i++ ) m_vals[i] = (*other)[i];
    }
    
    /// return the number of variables in this calibration
    inline unsigned int size() const {
      return m_vals.size();
    }
    
    /// set the status of the calibration
    inline void setStatus(STATUS stat) { m_status = stat; };
    /// return the status of the calibration
    inline STATUS getStatus() const { return m_status; }
    
    /// fill from a set of values
    void setVals(const std::vector<float>& vals, STATUS status);

    // fill from individaul values
    void setVals(float v1, float v2, STATUS stat);
    void setVals(float v1, float v2, float v3, STATUS stat);
    void setVals(float v1, float v2, float v3, float v4, STATUS stat);
    void setVals(float v1, float v2, float v3, float v4, float v5, STATUS stat);
    void setVals(float v1, float v2, float v3, float v4, float v5, float v6, STATUS stat);
    void setVals(float v1, float v2, float v3, float v4, float v5, float v6, float v7, STATUS stat);
    
  protected:
    
    /// I really don't know what this is doing here
    AcdCalibObj();

    /// set the calibration status
    bool setStatusInt(int stat);
    
    /// allocate space to store the values
    inline void sizeVals(int n) {
      m_vals.resize(n);
    }
    
  private:
    
    /// Status of the calibration
    STATUS m_status;
    
    /// Calibration values
    std::vector<float> m_vals;
  };
}


#endif

