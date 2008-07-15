// $Header$
#ifndef CalibData_AcdHighPed_h
#define CalibData_AcdHighPed_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdHighPedestalFitDesc
   *
   * @brief Description of an ACD pedestal calibration.
   * 
   * This calibration consists of:
   *  - mean  = the pedestal in PHA counts
   *  - width = the width of the pedestal
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdHighPedestalFitDesc : public AcdCalibDescription {
  public:
    /// Get this description
    static const AcdHighPedestalFitDesc& instance() {
      static const AcdHighPedestalFitDesc desc;
      return desc;
    }
  public:
    /// Trivial D'ctor
    virtual ~AcdHighPedestalFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdHighPedestalFitDesc()
      :AcdCalibDescription(AcdCalibData::PED_HIGH,"ACD_HighPed"){
      addVarName("mean");
      addVarName("width");
    }
  };

  /** 
   * @class AcdHighPed
   *
   * @brief An ACD pedestal calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - mean  = the pedestal in PHA counts
   *  - width = the width of the pedestal
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdHighPed : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_HighPed;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::PED_HIGH;
    }
  public:
    /// Build from description and a set of values
    AcdHighPed(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdHighPed(float mean, float width, STATUS status) :
      AcdCalibObj(status,AcdHighPedestalFitDesc::instance()){
      setVals(mean,width,status);
    }
    /// Trivial d'tor
    virtual ~AcdHighPed() {;}

    // Provide access to the values
    float getMean() const {return (*this)[0];}
    float getWidth() const {return (*this)[1]; }

  };
}

#endif
