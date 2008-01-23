// $Header$
#ifndef CalibData_AcdHighRange_h
#define CalibData_AcdHighRange_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdHighRangeFitDesc
   *
   * @brief Description of an ACD calibration for high range readout
   * 
   * This calibration consists of:
   *  - pedestal   = High range pedestal in PHA counts
   *  - slope      = Mips / PHA count near pedestal
   *  - saturation = Electronics saturation point in PHA counts
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdHighRangeFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdHighRangeFitDesc& instance() {
      static const AcdHighRangeFitDesc desc;
      return desc;
    }       
  public:
    /// Trivial D'ctor
    virtual ~AcdHighRangeFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdHighRangeFitDesc()
      :AcdCalibDescription(AcdCalibData::HIGH_RANGE,"ACD_HighRange"){
      addVarName("pedestal");
      addVarName("slope");
      addVarName("saturation");      
    }
  };

  /** 
   * @class AcdHighRange
   *
   * @brief An ACD calibration for high range readout for 1 PMT.
   * 
   * This calibration consists of:
   *  - pedestal   = High range pedestal in PHA counts
   *  - slope      = Mips / PHA count near pedestal
   *  - saturation = Electronics saturation point in PHA counts
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdHighRange : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_HighRange;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::HIGH_RANGE;
    }
  public:
    /// Build from description and a set of values
    AcdHighRange(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdHighRange(float pedestal, float slope, float saturation, STATUS status) :
      AcdCalibObj(status,AcdHighRangeFitDesc::instance()){
      setVals(pedestal,slope,saturation,status);
    }
    /// Trivial d'tor
    virtual ~AcdHighRange() {}

    // Provide access to the values
    float getPedestal() const { return (*this)[0];}
    float getSlope() const { return (*this)[1]; }
    float getSaturation() const { return (*this)[2]; }

  };
}


#endif
