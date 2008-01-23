// $Header$
#ifndef CalibData_AcdPed_h
#define CalibData_AcdPed_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdPedestalFitDesc
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

  class AcdPedestalFitDesc : public AcdCalibDescription {
  public:
    /// Get this description
    static const AcdPedestalFitDesc& instance() {
      static const AcdPedestalFitDesc desc;
      return desc;
    }
  public:
    /// Trivial D'ctor
    virtual ~AcdPedestalFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdPedestalFitDesc()
      :AcdCalibDescription(AcdCalibData::PEDESTAL,"ACD_Ped"){
      addVarName("mean");
      addVarName("width");
    }
  };

  /** 
   * @class AcdPed
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

  class AcdPed : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_Ped;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::PEDESTAL;
    }
  public:
    /// Build from description and a set of values
    AcdPed(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdPed(float mean, float width, STATUS status) :
      AcdCalibObj(status,AcdPedestalFitDesc::instance()){
      setVals(mean,width,status);
    }
    /// Trivial d'tor
    virtual ~AcdPed() {;}

    // Provide access to the values
    float getMean() const {return (*this)[0];}
    float getWidth() const {return (*this)[1]; }

  };
}

#endif
