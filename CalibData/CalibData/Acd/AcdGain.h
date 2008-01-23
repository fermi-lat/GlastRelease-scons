// $Header$
#ifndef CalibData_AcdGain_h
#define CalibData_AcdGain_h

#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibEnum.h"

namespace CalibData {

  /** 
   * @class AcdGainFitDesc
   *
   * @brief Description of Gain (aka MIP peak) calibration
   * 
   * This calibration consists of:
   *  - peak  = the mip peak in PHA counts above pedestal
   *  - width = the width of the MIP peak in PHA counts
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdGainFitDesc : public AcdCalibDescription {
  public:
    /// Get this description
   static const AcdGainFitDesc & instance() {
      static const AcdGainFitDesc desc;
      return desc;
    }
  public:
    /// Trivial D'ctor
    virtual ~AcdGainFitDesc(){;};    
  private:    
    /// This is a singleton
    AcdGainFitDesc()
      :AcdCalibDescription(AcdCalibData::GAIN,"ACD_Gain"){
      addVarName("peak");
      addVarName("width");
    }
  };

  /** 
   * @class AcdGain
   *
   * @brief A Gain (aka MIP peak) calibration for 1 PMT.
   * 
   * This calibration consists of:
   *  - peak  = the mip peak in PHA counts above pedestal
   *  - width = the width of the MIP peak in PHA counts
   *
   * @author Eric Charles
   * $Header$
   */

  class AcdGain : public AcdCalibObj {    
  public:
    /// For gaudi
    static const CLID& calibCLID() {
      return CLID_Calib_ACD_ElecGain;
    }
    /// Define the type of calibration
    static AcdCalibData::CALTYPE calibType() {
      return AcdCalibData::GAIN;
    }
  public:
    /// Build from description and a set of values
    AcdGain(const AcdCalibDescription& desc, const std::vector<float>& vals, STATUS status=NOFIT) :
      AcdCalibObj(status,vals,desc){
      assert( desc.calibType() == calibType() );
      setVals(vals,status);
    }
    /// Build from individaul values
    AcdGain(float peak, float width, STATUS status) :
      AcdCalibObj(status,AcdGainFitDesc::instance()){
      setVals(peak,width,status);
    }
     /// Trivial d'tor
    virtual ~AcdGain() {}

    // Provide access to the values
    float getPeak() const { return (*this)[0];}
    float getWidth() const { return (*this)[1]; }
  };
}

#endif
