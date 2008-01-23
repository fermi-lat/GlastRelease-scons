// $Header$
#ifndef CalibData_AcdRange_h
#define CalibData_AcdRange_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdRangeFitDesc
   *
   * @brief Description of an ACD range crossover calibration.
   * 
   * This calibration consists of:
   *  - low_max  = Highest PHA value read in LOW range
   *  - high_min = Lowest PHA value read in HIGH range, slightly above HIGH range pedestal
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdRangeFitDesc : public AcdCalibDescription {
  public:    
    /// Get this description
    static const AcdRangeFitDesc& instance() {
      static const AcdRangeFitDesc desc;
      return desc;
    };        
  public:
    /// Trivial D'ctor
    virtual ~AcdRangeFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdRangeFitDesc()
      :AcdCalibDescription(AcdCalibData::RANGE,"ACD_Range"){
      addVarName("low_max");
      addVarName("high_min");
    }
  };

 /** 
   * @class AcdRange
   *
   * @brief An ACD range crossover calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - low_max  = Highest PHA value read in LOW range
   *  - high_min = Lowest PHA value read in HIGH range, slightly above HIGH range pedestal
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdRange : public AcdCalibObj {
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_Range;
    } 
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::RANGE;
    }
  public:
    /// Build from description and a set of values
    AcdRange(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdRange(float lowMax, float highMin, STATUS status) :
      AcdCalibObj(status,AcdRangeFitDesc::instance()){
      setVals(lowMax,highMin,status);
    }
    /// Trivial d'tor
    virtual ~AcdRange() {}

    // Provide access to the values
    float getLowMax() const { return (*this)[0];}
    float getHighMin() const { return (*this)[1]; }
  };
}


#endif
