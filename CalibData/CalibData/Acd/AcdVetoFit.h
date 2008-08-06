// $Header$
#ifndef CalibData_AcdVetoFit_h
#define CalibData_AcdVetoFit_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdVetoFitFitDesc
   *
   * @brief Description of a VetoFit threshold calibration.
   * 
   * This calibration consists of:
   *  - veto  = the VetoFit threshold in PHA counts (50% point of turn on curve)
   *  - width = the width of the turn on curve
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdVetoFitFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdVetoFitFitDesc& instance() {
      static const AcdVetoFitFitDesc desc;
      return desc;
    };     
  public:
    /// Trivial D'ctor
    virtual ~AcdVetoFitFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdVetoFitFitDesc()
      :AcdCalibDescription(AcdCalibData::VETO_FIT,"ACD_VetoFit"){
      addVarName("slope");
      addVarName("offset");
      addVarName("mipPeak");
    }
  };

  /** 
   * @class AcdVetoFit
   *
   * @brief A VetoFit threshold calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - veto  = the VetoFit threshold in PHA counts (50% point of turn on curve)
   *  - width = the width of the turn on curve
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdVetoFit : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_VetoFit;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::VETO_FIT;
    }
  public:
    /// Build from description and a set of values
    AcdVetoFit(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdVetoFit(float veto, float width, STATUS status) :
      AcdCalibObj(status,AcdVetoFitFitDesc::instance()){
      setVals(veto,width,status);
    }
    /// Trivial d'tor
    virtual ~AcdVetoFit() {}

    // Provide access to the values
    float getSlope() const { return (*this)[0];}
    float getOffset() const { return (*this)[1]; }
    float getMipPeak() const { return (*this)[2]; }
  };
}


#endif
